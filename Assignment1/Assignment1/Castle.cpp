//***************************************************************************************
// LitColumnsApp.cpp by Frank Luna (C) 2015 All Rights Reserved.
//***************************************************************************************

#include "../../Common/d3dApp.h"
#include "../../Common/MathHelper.h"
#include "../../Common/UploadBuffer.h"
#include "GeometryGenerator.h"
#include "FrameResource.h"

using Microsoft::WRL::ComPtr;
using namespace DirectX;
using namespace DirectX::PackedVector;

#pragma comment(lib, "d3dcompiler.lib")
#pragma comment(lib, "D3D12.lib")

const int gNumFrameResources = 3;

// Lightweight structure stores parameters to draw a shape.  This will
// vary from app-to-app.
struct RenderItem
{
	RenderItem() = default;

	// World matrix of the shape that describes the object's local space
	// relative to the world space, which defines the position, orientation,
	// and scale of the object in the world.
	XMFLOAT4X4 World = MathHelper::Identity4x4();

	XMFLOAT4X4 TexTransform = MathHelper::Identity4x4();

	// Dirty flag indicating the object data has changed and we need to update the constant buffer.
	// Because we have an object cbuffer for each FrameResource, we have to apply the
	// update to each FrameResource.  Thus, when we modify obect data we should set 
	// NumFramesDirty = gNumFrameResources so that each frame resource gets the update.
	int NumFramesDirty = gNumFrameResources;

	// Index into GPU constant buffer corresponding to the ObjectCB for this render item.
	UINT ObjCBIndex = -1;

	Material* Mat = nullptr;
	MeshGeometry* Geo = nullptr;

	// Primitive topology.
	D3D12_PRIMITIVE_TOPOLOGY PrimitiveType = D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST;

	// DrawIndexedInstanced parameters.
	UINT IndexCount = 0;
	UINT StartIndexLocation = 0;
	int BaseVertexLocation = 0;
};

class LitColumnsApp : public D3DApp
{
public:
	LitColumnsApp(HINSTANCE hInstance);
	LitColumnsApp(const LitColumnsApp& rhs) = delete;
	LitColumnsApp& operator=(const LitColumnsApp& rhs) = delete;
	~LitColumnsApp();

	virtual bool Initialize()override;

private:
	virtual void OnResize()override;
	virtual void Update(const GameTimer& gt)override;
	virtual void Draw(const GameTimer& gt)override;

	virtual void OnMouseDown(WPARAM btnState, int x, int y)override;
	virtual void OnMouseUp(WPARAM btnState, int x, int y)override;
	virtual void OnMouseMove(WPARAM btnState, int x, int y)override;

	void OnKeyboardInput(const GameTimer& gt);
	void UpdateCamera(const GameTimer& gt);
	void AnimateMaterials(const GameTimer& gt);
	void UpdateObjectCBs(const GameTimer& gt);
	void UpdateMaterialCBs(const GameTimer& gt);
	void UpdateMainPassCB(const GameTimer& gt);

	void BuildRootSignature();
	void BuildShadersAndInputLayout();
	void BuildShapeGeometry();
	void BuildPSOs();
	void BuildFrameResources();
	void BuildMaterials();
	void BuildRenderItem(int index, std::string itemName, std::string material, XMFLOAT3 scaling, XMFLOAT3 rotation, XMFLOAT3 translation);
	void BuildRenderItems();
	void DrawRenderItems(ID3D12GraphicsCommandList* cmdList, const std::vector<RenderItem*>& ritems);

private:

	std::vector<std::unique_ptr<FrameResource>> mFrameResources;
	FrameResource* mCurrFrameResource = nullptr;
	int mCurrFrameResourceIndex = 0;

	UINT mCbvSrvDescriptorSize = 0;

	ComPtr<ID3D12RootSignature> mRootSignature = nullptr;

	ComPtr<ID3D12DescriptorHeap> mSrvDescriptorHeap = nullptr;

	std::unordered_map<std::string, std::unique_ptr<MeshGeometry>> mGeometries;
	std::unordered_map<std::string, std::unique_ptr<Material>> mMaterials;
	std::unordered_map<std::string, std::unique_ptr<Texture>> mTextures;
	std::unordered_map<std::string, ComPtr<ID3DBlob>> mShaders;

	std::vector<D3D12_INPUT_ELEMENT_DESC> mInputLayout;

	ComPtr<ID3D12PipelineState> mOpaquePSO = nullptr;

	// List of all the render items.
	std::vector<std::unique_ptr<RenderItem>> mAllRitems;

	// Render items divided by PSO.
	std::vector<RenderItem*> mOpaqueRitems;

	PassConstants mMainPassCB;

	XMFLOAT3 mEyePos = { 85.0f,  51.0f,  40.0f };
	XMFLOAT4X4 mView = MathHelper::Identity4x4();
	XMFLOAT4X4 mProj = MathHelper::Identity4x4();

	float mTheta = 1.5f*XM_PI;
	float mPhi = 0.2f*XM_PI;
	float mRadius = 15.0f;

	POINT mLastMousePos;
};

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE prevInstance,
	PSTR cmdLine, int showCmd)
{
	// Enable run-time memory check for debug builds.
#if defined(DEBUG) | defined(_DEBUG)
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif

	try
	{
		LitColumnsApp theApp(hInstance);
		if (!theApp.Initialize())
			return 0;

		return theApp.Run();
	}
	catch (DxException& e)
	{
		MessageBox(nullptr, e.ToString().c_str(), L"HR Failed", MB_OK);
		return 0;
	}
}

LitColumnsApp::LitColumnsApp(HINSTANCE hInstance)
	: D3DApp(hInstance)
{
}

LitColumnsApp::~LitColumnsApp()
{
	if (md3dDevice != nullptr)
		FlushCommandQueue();
}

bool LitColumnsApp::Initialize()
{
	if (!D3DApp::Initialize())
		return false;

	// Reset the command list to prep for initialization commands.
	ThrowIfFailed(mCommandList->Reset(mDirectCmdListAlloc.Get(), nullptr));

	// Get the increment size of a descriptor in this heap type.  This is hardware specific, 
	// so we have to query this information.
	mCbvSrvDescriptorSize = md3dDevice->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV);

	BuildRootSignature();
	BuildShadersAndInputLayout();
	BuildShapeGeometry();
	BuildMaterials();
	BuildRenderItems();
	BuildFrameResources();
	BuildPSOs();

	// Execute the initialization commands.
	ThrowIfFailed(mCommandList->Close());
	ID3D12CommandList* cmdsLists[] = { mCommandList.Get() };
	mCommandQueue->ExecuteCommandLists(_countof(cmdsLists), cmdsLists);

	// Wait until initialization is complete.
	FlushCommandQueue();

	return true;
}

void LitColumnsApp::OnResize()
{
	D3DApp::OnResize();

	// The window resized, so update the aspect ratio and recompute the projection matrix.
	XMMATRIX P = XMMatrixPerspectiveFovLH(0.25f*MathHelper::Pi, AspectRatio(), 1.0f, 1000.0f);
	XMStoreFloat4x4(&mProj, P);
}

void LitColumnsApp::Update(const GameTimer& gt)
{
	OnKeyboardInput(gt);
	UpdateCamera(gt);

	// Cycle through the circular frame resource array.
	mCurrFrameResourceIndex = (mCurrFrameResourceIndex + 1) % gNumFrameResources;
	mCurrFrameResource = mFrameResources[mCurrFrameResourceIndex].get();

	// Has the GPU finished processing the commands of the current frame resource?
	// If not, wait until the GPU has completed commands up to this fence point.
	if (mCurrFrameResource->Fence != 0 && mFence->GetCompletedValue() < mCurrFrameResource->Fence)
	{
		HANDLE eventHandle = CreateEventEx(nullptr, false, false, EVENT_ALL_ACCESS);
		ThrowIfFailed(mFence->SetEventOnCompletion(mCurrFrameResource->Fence, eventHandle));
		WaitForSingleObject(eventHandle, INFINITE);
		CloseHandle(eventHandle);
	}

	AnimateMaterials(gt);
	UpdateObjectCBs(gt);
	UpdateMaterialCBs(gt);
	UpdateMainPassCB(gt);
}

void LitColumnsApp::Draw(const GameTimer& gt)
{
	auto cmdListAlloc = mCurrFrameResource->CmdListAlloc;

	// Reuse the memory associated with command recording.
	// We can only reset when the associated command lists have finished execution on the GPU.
	ThrowIfFailed(cmdListAlloc->Reset());

	// A command list can be reset after it has been added to the command queue via ExecuteCommandList.
	// Reusing the command list reuses memory.
	ThrowIfFailed(mCommandList->Reset(cmdListAlloc.Get(), mOpaquePSO.Get()));

	mCommandList->RSSetViewports(1, &mScreenViewport);
	mCommandList->RSSetScissorRects(1, &mScissorRect);

	// Indicate a state transition on the resource usage.
	mCommandList->ResourceBarrier(1, &CD3DX12_RESOURCE_BARRIER::Transition(CurrentBackBuffer(),
		D3D12_RESOURCE_STATE_PRESENT, D3D12_RESOURCE_STATE_RENDER_TARGET));

	// Clear the back buffer and depth buffer.
	mCommandList->ClearRenderTargetView(CurrentBackBufferView(), Colors::LightSteelBlue, 0, nullptr);
	mCommandList->ClearDepthStencilView(DepthStencilView(), D3D12_CLEAR_FLAG_DEPTH | D3D12_CLEAR_FLAG_STENCIL, 1.0f, 0, 0, nullptr);

	// Specify the buffers we are going to render to.
	mCommandList->OMSetRenderTargets(1, &CurrentBackBufferView(), true, &DepthStencilView());

	mCommandList->SetGraphicsRootSignature(mRootSignature.Get());

	auto passCB = mCurrFrameResource->PassCB->Resource();
	mCommandList->SetGraphicsRootConstantBufferView(2, passCB->GetGPUVirtualAddress());

	DrawRenderItems(mCommandList.Get(), mOpaqueRitems);

	// Indicate a state transition on the resource usage.
	mCommandList->ResourceBarrier(1, &CD3DX12_RESOURCE_BARRIER::Transition(CurrentBackBuffer(),
		D3D12_RESOURCE_STATE_RENDER_TARGET, D3D12_RESOURCE_STATE_PRESENT));

	// Done recording commands.
	ThrowIfFailed(mCommandList->Close());

	// Add the command list to the queue for execution.
	ID3D12CommandList* cmdsLists[] = { mCommandList.Get() };
	mCommandQueue->ExecuteCommandLists(_countof(cmdsLists), cmdsLists);

	// Swap the back and front buffers
	ThrowIfFailed(mSwapChain->Present(0, 0));
	mCurrBackBuffer = (mCurrBackBuffer + 1) % SwapChainBufferCount;

	// Advance the fence value to mark commands up to this fence point.
	mCurrFrameResource->Fence = ++mCurrentFence;

	// Add an instruction to the command queue to set a new fence point. 
	// Because we are on the GPU timeline, the new fence point won't be 
	// set until the GPU finishes processing all the commands prior to this Signal().
	mCommandQueue->Signal(mFence.Get(), mCurrentFence);
}

void LitColumnsApp::OnMouseDown(WPARAM btnState, int x, int y)
{
	mLastMousePos.x = x;
	mLastMousePos.y = y;

	SetCapture(mhMainWnd);
}

void LitColumnsApp::OnMouseUp(WPARAM btnState, int x, int y)
{
	ReleaseCapture();
}

void LitColumnsApp::OnMouseMove(WPARAM btnState, int x, int y)
{
	if ((btnState & MK_LBUTTON) != 0)
	{
		// Make each pixel correspond to a quarter of a degree.
		float dx = XMConvertToRadians(0.25f*static_cast<float>(x - mLastMousePos.x));
		float dy = XMConvertToRadians(0.25f*static_cast<float>(y - mLastMousePos.y));

		// Update angles based on input to orbit camera around box.
		mTheta += dx;
		mPhi += dy;

		// Restrict the angle mPhi.
		mPhi = MathHelper::Clamp(mPhi, 0.1f, MathHelper::Pi - 0.1f);
	}
	else if ((btnState & MK_RBUTTON) != 0)
	{
		// Make each pixel correspond to 0.2 unit in the scene.
		float dx = 0.05f*static_cast<float>(x - mLastMousePos.x);
		float dy = 0.05f*static_cast<float>(y - mLastMousePos.y);

		// Update the camera radius based on input.
		mRadius += dx - dy;

		// Restrict the radius.
		mRadius = MathHelper::Clamp(mRadius, 5.0f, 150.0f);
	}

	mLastMousePos.x = x;
	mLastMousePos.y = y;
}

void LitColumnsApp::OnKeyboardInput(const GameTimer& gt)
{
}

void LitColumnsApp::UpdateCamera(const GameTimer& gt)
{
	// Convert Spherical to Cartesian coordinates.
	mEyePos.x = mRadius * sinf(mPhi)*cosf(mTheta);
	mEyePos.z = mRadius * sinf(mPhi)*sinf(mTheta);
	mEyePos.y = mRadius * cosf(mPhi);

	// Build the view matrix.
	XMVECTOR pos = XMVectorSet(mEyePos.x, mEyePos.y, mEyePos.z, 1.0f);
	XMVECTOR target = XMVectorZero();
	XMVECTOR up = XMVectorSet(0.0f, 1.0f, 0.0f, 0.0f);

	XMMATRIX view = XMMatrixLookAtLH(pos, target, up);
	XMStoreFloat4x4(&mView, view);
}

void LitColumnsApp::AnimateMaterials(const GameTimer& gt)
{

}

void LitColumnsApp::UpdateObjectCBs(const GameTimer& gt)
{
	auto currObjectCB = mCurrFrameResource->ObjectCB.get();
	for (auto& e : mAllRitems)
	{
		// Only update the cbuffer data if the constants have changed.  
		// This needs to be tracked per frame resource.
		if (e->NumFramesDirty > 0)
		{
			XMMATRIX world = XMLoadFloat4x4(&e->World);
			XMMATRIX texTransform = XMLoadFloat4x4(&e->TexTransform);

			ObjectConstants objConstants;
			XMStoreFloat4x4(&objConstants.World, XMMatrixTranspose(world));
			XMStoreFloat4x4(&objConstants.TexTransform, XMMatrixTranspose(texTransform));

			currObjectCB->CopyData(e->ObjCBIndex, objConstants);

			// Next FrameResource need to be updated too.
			e->NumFramesDirty--;
		}
	}
}

void LitColumnsApp::UpdateMaterialCBs(const GameTimer& gt)
{
	auto currMaterialCB = mCurrFrameResource->MaterialCB.get();
	for (auto& e : mMaterials)
	{
		// Only update the cbuffer data if the constants have changed.  If the cbuffer
		// data changes, it needs to be updated for each FrameResource.
		Material* mat = e.second.get();
		if (mat->NumFramesDirty > 0)
		{
			XMMATRIX matTransform = XMLoadFloat4x4(&mat->MatTransform);

			MaterialConstants matConstants;
			matConstants.DiffuseAlbedo = mat->DiffuseAlbedo;
			matConstants.FresnelR0 = mat->FresnelR0;
			matConstants.Roughness = mat->Roughness;
			XMStoreFloat4x4(&matConstants.MatTransform, XMMatrixTranspose(matTransform));

			currMaterialCB->CopyData(mat->MatCBIndex, matConstants);

			// Next FrameResource need to be updated too.
			mat->NumFramesDirty--;
		}
	}
}

void LitColumnsApp::UpdateMainPassCB(const GameTimer& gt)
{
	XMMATRIX view = XMLoadFloat4x4(&mView);
	XMMATRIX proj = XMLoadFloat4x4(&mProj);

	XMMATRIX viewProj = XMMatrixMultiply(view, proj);
	XMMATRIX invView = XMMatrixInverse(&XMMatrixDeterminant(view), view);
	XMMATRIX invProj = XMMatrixInverse(&XMMatrixDeterminant(proj), proj);
	XMMATRIX invViewProj = XMMatrixInverse(&XMMatrixDeterminant(viewProj), viewProj);

	XMStoreFloat4x4(&mMainPassCB.View, XMMatrixTranspose(view));
	XMStoreFloat4x4(&mMainPassCB.InvView, XMMatrixTranspose(invView));
	XMStoreFloat4x4(&mMainPassCB.Proj, XMMatrixTranspose(proj));
	XMStoreFloat4x4(&mMainPassCB.InvProj, XMMatrixTranspose(invProj));
	XMStoreFloat4x4(&mMainPassCB.ViewProj, XMMatrixTranspose(viewProj));
	XMStoreFloat4x4(&mMainPassCB.InvViewProj, XMMatrixTranspose(invViewProj));
	mMainPassCB.EyePosW = mEyePos;
	mMainPassCB.RenderTargetSize = XMFLOAT2((float)mClientWidth, (float)mClientHeight);
	mMainPassCB.InvRenderTargetSize = XMFLOAT2(1.0f / mClientWidth, 1.0f / mClientHeight);
	mMainPassCB.NearZ = 1.0f;
	mMainPassCB.FarZ = 1000.0f;
	mMainPassCB.TotalTime = gt.TotalTime();
	mMainPassCB.DeltaTime = gt.DeltaTime();
	mMainPassCB.AmbientLight = { 0.8f, 0.8f, 0.8f, 1.0f };
	auto currPassCB = mCurrFrameResource->PassCB.get();
	currPassCB->CopyData(0, mMainPassCB);
}

void LitColumnsApp::BuildRootSignature()
{
	// Root parameter can be a table, root descriptor or root constants.
	CD3DX12_ROOT_PARAMETER slotRootParameter[3];

	// Create root CBV.
	slotRootParameter[0].InitAsConstantBufferView(0);
	slotRootParameter[1].InitAsConstantBufferView(1);
	slotRootParameter[2].InitAsConstantBufferView(2);

	// A root signature is an array of root parameters.
	CD3DX12_ROOT_SIGNATURE_DESC rootSigDesc(3, slotRootParameter, 0, nullptr, D3D12_ROOT_SIGNATURE_FLAG_ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT);

	// create a root signature with a single slot which points to a descriptor range consisting of a single constant buffer
	ComPtr<ID3DBlob> serializedRootSig = nullptr;
	ComPtr<ID3DBlob> errorBlob = nullptr;
	HRESULT hr = D3D12SerializeRootSignature(&rootSigDesc, D3D_ROOT_SIGNATURE_VERSION_1,
		serializedRootSig.GetAddressOf(), errorBlob.GetAddressOf());

	if (errorBlob != nullptr)
	{
		::OutputDebugStringA((char*)errorBlob->GetBufferPointer());
	}
	ThrowIfFailed(hr);

	ThrowIfFailed(md3dDevice->CreateRootSignature(
		0,
		serializedRootSig->GetBufferPointer(),
		serializedRootSig->GetBufferSize(),
		IID_PPV_ARGS(mRootSignature.GetAddressOf())));
}

void LitColumnsApp::BuildShadersAndInputLayout()
{
	const D3D_SHADER_MACRO alphaTestDefines[] =
	{
		"ALPHA_TEST", "1",
		NULL, NULL
	};

	mShaders["standardVS"] = d3dUtil::CompileShader(L"Shaders\\Default.hlsl", nullptr, "VS", "vs_5_1");
	mShaders["opaquePS"] = d3dUtil::CompileShader(L"Shaders\\Default.hlsl", nullptr, "PS", "ps_5_1");

	mInputLayout =
	{
		{ "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
		{ "NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
		{ "TEXCOORD", 0, DXGI_FORMAT_R32G32_FLOAT, 0, 24, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
	};
}

void LitColumnsApp::BuildShapeGeometry()
{
	GeometryGenerator geoGen;
	GeometryGenerator::MeshData box = geoGen.CreateBox(1.5f, 1.5f, 1.5f, 3);
	GeometryGenerator::MeshData grid = geoGen.CreateGrid(20.0f, 30.0f, 60, 40);
	GeometryGenerator::MeshData sphere = geoGen.CreateSphere(0.5f, 20, 20);
	GeometryGenerator::MeshData cylinder = geoGen.CreateCylinder(0.5f, 0.0f, 3.0f, 20, 20);

	//// it all starts now
	GeometryGenerator::MeshData star = geoGen.CreateStar();
	GeometryGenerator::MeshData wedge = geoGen.CreateWedge();
	GeometryGenerator::MeshData triPrism = geoGen.CreateTriPrism();
	GeometryGenerator::MeshData cube = geoGen.CreateCube();
	GeometryGenerator::MeshData pyramid = geoGen.CreatePyramid();
	GeometryGenerator::MeshData hexagon = geoGen.CreateHexagon();
	GeometryGenerator::MeshData simsInd = geoGen.CreateSimsIndicator();
	GeometryGenerator::MeshData myCylinder = geoGen.CreateMyCylinder();
	GeometryGenerator::MeshData cone = geoGen.CreateCone();
	GeometryGenerator::MeshData mySphere = geoGen.CreateMySphere();

	//
	// We are concatenating all the geometry into one big vertex/index buffer.  So
	// define the regions in the buffer each submesh covers.
	//

	// Cache the vertex offsets to each object in the concatenated vertex buffer.
	UINT boxVertexOffset = 0;
	UINT gridVertexOffset = (UINT)box.Vertices.size();
	UINT sphereVertexOffset = gridVertexOffset + (UINT)grid.Vertices.size();
	UINT cylinderVertexOffset = sphereVertexOffset + (UINT)sphere.Vertices.size();

	UINT starVertexOffset = cylinderVertexOffset + (UINT)cylinder.Vertices.size();
	UINT wedgeVertexOffset = starVertexOffset + (UINT)star.Vertices.size();
	UINT triPrismVertexOffset = wedgeVertexOffset + (UINT)wedge.Vertices.size();
	UINT cubeVertexOffset = triPrismVertexOffset + (UINT)triPrism.Vertices.size();
	UINT pyramidVertexOffset = cubeVertexOffset + (UINT)cube.Vertices.size();
	UINT hexagonVertexOffset = pyramidVertexOffset + (UINT)pyramid.Vertices.size();
	UINT simsIndVertexOffset = hexagonVertexOffset + (UINT)hexagon.Vertices.size();
	UINT myCylinderVertexOffset = simsIndVertexOffset + (UINT)simsInd.Vertices.size();
	UINT coneVertexOffset = myCylinderVertexOffset + (UINT)myCylinder.Vertices.size();
	UINT mySphereVertexOffset = coneVertexOffset + (UINT)cone.Vertices.size();

	// Cache the starting index for each object in the concatenated index buffer.
	UINT boxIndexOffset = 0;
	UINT gridIndexOffset = (UINT)box.Indices32.size();
	UINT sphereIndexOffset = gridIndexOffset + (UINT)grid.Indices32.size();
	UINT cylinderIndexOffset = sphereIndexOffset + (UINT)sphere.Indices32.size();

	UINT starIndexOffset = cylinderIndexOffset + (UINT)cylinder.Indices32.size();
	UINT wedgeIndexOffset = starIndexOffset + (UINT)star.Indices32.size();
	UINT triPrismIndexOffset = wedgeIndexOffset + (UINT)wedge.Indices32.size();
	UINT cubeIndexOffset = triPrismIndexOffset + (UINT)triPrism.Indices32.size();
	UINT pyramidIndexOffset = cubeIndexOffset + (UINT)cube.Indices32.size();
	UINT hexagonIndexOffset = pyramidIndexOffset + (UINT)pyramid.Indices32.size();
	UINT simsIndIndexOffset = hexagonIndexOffset + (UINT)hexagon.Indices32.size();
	UINT myCylinderIndexOffset = simsIndIndexOffset + (UINT)simsInd.Indices32.size();
	UINT coneIndexOffset = myCylinderIndexOffset + (UINT)myCylinder.Indices32.size();
	UINT mySphereIndexOffset = coneIndexOffset + (UINT)cone.Indices32.size();

	SubmeshGeometry boxSubmesh;
	boxSubmesh.IndexCount = (UINT)box.Indices32.size();
	boxSubmesh.StartIndexLocation = boxIndexOffset;
	boxSubmesh.BaseVertexLocation = boxVertexOffset;

	SubmeshGeometry gridSubmesh;
	gridSubmesh.IndexCount = (UINT)grid.Indices32.size();
	gridSubmesh.StartIndexLocation = gridIndexOffset;
	gridSubmesh.BaseVertexLocation = gridVertexOffset;

	SubmeshGeometry sphereSubmesh;
	sphereSubmesh.IndexCount = (UINT)sphere.Indices32.size();
	sphereSubmesh.StartIndexLocation = sphereIndexOffset;
	sphereSubmesh.BaseVertexLocation = sphereVertexOffset;

	SubmeshGeometry cylinderSubmesh;
	cylinderSubmesh.IndexCount = (UINT)cylinder.Indices32.size();
	cylinderSubmesh.StartIndexLocation = cylinderIndexOffset;
	cylinderSubmesh.BaseVertexLocation = cylinderVertexOffset;

	SubmeshGeometry starSubmesh;
	starSubmesh.IndexCount = (UINT)star.Indices32.size();
	starSubmesh.StartIndexLocation = starIndexOffset;
	starSubmesh.BaseVertexLocation = starVertexOffset;

	SubmeshGeometry wedgeSubmesh;
	wedgeSubmesh.IndexCount = (UINT)wedge.Indices32.size();
	wedgeSubmesh.StartIndexLocation = wedgeIndexOffset;
	wedgeSubmesh.BaseVertexLocation = wedgeVertexOffset;

	SubmeshGeometry triPrismSubmesh;
	triPrismSubmesh.IndexCount = (UINT)triPrism.Indices32.size();
	triPrismSubmesh.StartIndexLocation = triPrismIndexOffset;
	triPrismSubmesh.BaseVertexLocation = triPrismVertexOffset;

	SubmeshGeometry cubeSubmesh;
	cubeSubmesh.IndexCount = (UINT)cube.Indices32.size();
	cubeSubmesh.StartIndexLocation = cubeIndexOffset;
	cubeSubmesh.BaseVertexLocation = cubeVertexOffset;

	SubmeshGeometry pyramidSubmesh;
	pyramidSubmesh.IndexCount = (UINT)pyramid.Indices32.size();
	pyramidSubmesh.StartIndexLocation = pyramidIndexOffset;
	pyramidSubmesh.BaseVertexLocation = pyramidVertexOffset;

	SubmeshGeometry hexagonSubmesh;
	hexagonSubmesh.IndexCount = (UINT)hexagon.Indices32.size();
	hexagonSubmesh.StartIndexLocation = hexagonIndexOffset;
	hexagonSubmesh.BaseVertexLocation = hexagonVertexOffset;

	SubmeshGeometry simsIndSubmesh;
	simsIndSubmesh.IndexCount = (UINT)simsInd.Indices32.size();
	simsIndSubmesh.StartIndexLocation = simsIndIndexOffset;
	simsIndSubmesh.BaseVertexLocation = simsIndVertexOffset;

	SubmeshGeometry myCylinderSubmesh;
	myCylinderSubmesh.IndexCount = (UINT)myCylinder.Indices32.size();
	myCylinderSubmesh.StartIndexLocation = myCylinderIndexOffset;
	myCylinderSubmesh.BaseVertexLocation = myCylinderVertexOffset;

	SubmeshGeometry coneSubmesh;
	coneSubmesh.IndexCount = (UINT)cone.Indices32.size();
	coneSubmesh.StartIndexLocation = coneIndexOffset;
	coneSubmesh.BaseVertexLocation = coneVertexOffset;

	SubmeshGeometry mySphereSubmesh;
	mySphereSubmesh.IndexCount = (UINT)mySphere.Indices32.size();
	mySphereSubmesh.StartIndexLocation = mySphereIndexOffset;
	mySphereSubmesh.BaseVertexLocation = mySphereVertexOffset;
	//
	// Extract the vertex elements we are interested in and pack the
	// vertices of all the meshes into one vertex buffer.
	//

	auto totalVertexCount =
		box.Vertices.size() +
		grid.Vertices.size() +
		sphere.Vertices.size() +
		cylinder.Vertices.size() +
		star.Vertices.size() +
		wedge.Vertices.size() +
		triPrism.Vertices.size() +
		cube.Vertices.size() +
		pyramid.Vertices.size() +
		hexagon.Vertices.size() +
		simsInd.Vertices.size() +
		myCylinder.Vertices.size() +
		cone.Vertices.size() +
		mySphere.Vertices.size();

	std::vector<Vertex> vertices(totalVertexCount);

	UINT k = 0;
	for (size_t i = 0; i < box.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = box.Vertices[i].Position;
		vertices[k].Normal = box.Vertices[i].Normal;
	}

	for (size_t i = 0; i < grid.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = grid.Vertices[i].Position;
		vertices[k].Normal = grid.Vertices[i].Normal;
	}

	for (size_t i = 0; i < sphere.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = sphere.Vertices[i].Position;
		vertices[k].Normal = sphere.Vertices[i].Normal;
	}

	for (size_t i = 0; i < cylinder.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = cylinder.Vertices[i].Position;
		vertices[k].Normal = cylinder.Vertices[i].Normal;
	}


	for (size_t i = 0; i < star.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = star.Vertices[i].Position;
		vertices[k].Normal = star.Vertices[i].Normal;
	}

	for (size_t i = 0; i < wedge.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = wedge.Vertices[i].Position;
		vertices[k].Normal = wedge.Vertices[i].Normal;
	}

	for (size_t i = 0; i < triPrism.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = triPrism.Vertices[i].Position;
		vertices[k].Normal = triPrism.Vertices[i].Normal;
	}

	for (size_t i = 0; i < cube.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = cube.Vertices[i].Position;
		vertices[k].Normal = cube.Vertices[i].Normal;
	}

	for (size_t i = 0; i < pyramid.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = pyramid.Vertices[i].Position;
		vertices[k].Normal = pyramid.Vertices[i].Normal;
	}

	for (size_t i = 0; i < hexagon.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = hexagon.Vertices[i].Position;
		vertices[k].Normal = hexagon.Vertices[i].Normal;
	}

	for (size_t i = 0; i < simsInd.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = simsInd.Vertices[i].Position;
		vertices[k].Normal = simsInd.Vertices[i].Normal;
	}

	for (size_t i = 0; i < myCylinder.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = myCylinder.Vertices[i].Position;
		vertices[k].Normal = myCylinder.Vertices[i].Normal;
	}

	for (size_t i = 0; i < cone.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = cone.Vertices[i].Position;
		vertices[k].Normal = cone.Vertices[i].Normal;
	}

	for (size_t i = 0; i < mySphere.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = mySphere.Vertices[i].Position;
		vertices[k].Normal = mySphere.Vertices[i].Normal;
	}

	std::vector<std::uint16_t> indices;
	indices.insert(indices.end(), std::begin(box.GetIndices16()), std::end(box.GetIndices16()));
	indices.insert(indices.end(), std::begin(grid.GetIndices16()), std::end(grid.GetIndices16()));
	indices.insert(indices.end(), std::begin(sphere.GetIndices16()), std::end(sphere.GetIndices16()));
	indices.insert(indices.end(), std::begin(cylinder.GetIndices16()), std::end(cylinder.GetIndices16()));

	indices.insert(indices.end(), std::begin(star.GetIndices16()), std::end(star.GetIndices16()));
	indices.insert(indices.end(), std::begin(wedge.GetIndices16()), std::end(wedge.GetIndices16()));
	indices.insert(indices.end(), std::begin(triPrism.GetIndices16()), std::end(triPrism.GetIndices16()));
	indices.insert(indices.end(), std::begin(cube.GetIndices16()), std::end(cube.GetIndices16()));
	indices.insert(indices.end(), std::begin(pyramid.GetIndices16()), std::end(pyramid.GetIndices16()));
	indices.insert(indices.end(), std::begin(hexagon.GetIndices16()), std::end(hexagon.GetIndices16()));
	indices.insert(indices.end(), std::begin(simsInd.GetIndices16()), std::end(simsInd.GetIndices16()));
	indices.insert(indices.end(), std::begin(myCylinder.GetIndices16()), std::end(myCylinder.GetIndices16()));
	indices.insert(indices.end(), std::begin(cone.GetIndices16()), std::end(cone.GetIndices16()));
	indices.insert(indices.end(), std::begin(mySphere.GetIndices16()), std::end(mySphere.GetIndices16()));

	const UINT vbByteSize = (UINT)vertices.size() * sizeof(Vertex);
	const UINT ibByteSize = (UINT)indices.size() * sizeof(std::uint16_t);

	auto geo = std::make_unique<MeshGeometry>();
	geo->Name = "shapeGeo";

	ThrowIfFailed(D3DCreateBlob(vbByteSize, &geo->VertexBufferCPU));
	CopyMemory(geo->VertexBufferCPU->GetBufferPointer(), vertices.data(), vbByteSize);

	ThrowIfFailed(D3DCreateBlob(ibByteSize, &geo->IndexBufferCPU));
	CopyMemory(geo->IndexBufferCPU->GetBufferPointer(), indices.data(), ibByteSize);

	geo->VertexBufferGPU = d3dUtil::CreateDefaultBuffer(md3dDevice.Get(),
		mCommandList.Get(), vertices.data(), vbByteSize, geo->VertexBufferUploader);

	geo->IndexBufferGPU = d3dUtil::CreateDefaultBuffer(md3dDevice.Get(),
		mCommandList.Get(), indices.data(), ibByteSize, geo->IndexBufferUploader);

	geo->VertexByteStride = sizeof(Vertex);
	geo->VertexBufferByteSize = vbByteSize;
	geo->IndexFormat = DXGI_FORMAT_R16_UINT;
	geo->IndexBufferByteSize = ibByteSize;

	geo->DrawArgs["box"] = boxSubmesh;
	geo->DrawArgs["grid"] = gridSubmesh;
	geo->DrawArgs["sphere"] = sphereSubmesh;
	geo->DrawArgs["cylinder"] = cylinderSubmesh;

	geo->DrawArgs["star"] = starSubmesh;
	geo->DrawArgs["wedge"] = wedgeSubmesh;
	geo->DrawArgs["triPrism"] = triPrismSubmesh;
	geo->DrawArgs["cube"] = cubeSubmesh;
	geo->DrawArgs["pyramid"] = pyramidSubmesh;
	geo->DrawArgs["hexagon"] = hexagonSubmesh;
	geo->DrawArgs["simsInd"] = simsIndSubmesh;
	geo->DrawArgs["myCylinder"] = myCylinderSubmesh;
	geo->DrawArgs["cone"] = coneSubmesh;
	geo->DrawArgs["mySphere"] = mySphereSubmesh;

	mGeometries[geo->Name] = std::move(geo);
}

void LitColumnsApp::BuildPSOs()
{
	D3D12_GRAPHICS_PIPELINE_STATE_DESC opaquePsoDesc;

	//
	// PSO for opaque objects.
	//
	ZeroMemory(&opaquePsoDesc, sizeof(D3D12_GRAPHICS_PIPELINE_STATE_DESC));
	opaquePsoDesc.InputLayout = { mInputLayout.data(), (UINT)mInputLayout.size() };
	opaquePsoDesc.pRootSignature = mRootSignature.Get();
	opaquePsoDesc.VS =
	{
		reinterpret_cast<BYTE*>(mShaders["standardVS"]->GetBufferPointer()),
		mShaders["standardVS"]->GetBufferSize()
	};
	opaquePsoDesc.PS =
	{
		reinterpret_cast<BYTE*>(mShaders["opaquePS"]->GetBufferPointer()),
		mShaders["opaquePS"]->GetBufferSize()
	};
	opaquePsoDesc.RasterizerState = CD3DX12_RASTERIZER_DESC(D3D12_DEFAULT);
	opaquePsoDesc.BlendState = CD3DX12_BLEND_DESC(D3D12_DEFAULT);
	opaquePsoDesc.DepthStencilState = CD3DX12_DEPTH_STENCIL_DESC(D3D12_DEFAULT);
	opaquePsoDesc.SampleMask = UINT_MAX;
	opaquePsoDesc.PrimitiveTopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;
	opaquePsoDesc.NumRenderTargets = 1;
	opaquePsoDesc.RTVFormats[0] = mBackBufferFormat;
	opaquePsoDesc.SampleDesc.Count = m4xMsaaState ? 4 : 1;
	opaquePsoDesc.SampleDesc.Quality = m4xMsaaState ? (m4xMsaaQuality - 1) : 0;
	opaquePsoDesc.DSVFormat = mDepthStencilFormat;
	ThrowIfFailed(md3dDevice->CreateGraphicsPipelineState(&opaquePsoDesc, IID_PPV_ARGS(&mOpaquePSO)));
}

void LitColumnsApp::BuildFrameResources()
{
	for (int i = 0; i < gNumFrameResources; ++i)
	{
		mFrameResources.push_back(std::make_unique<FrameResource>(md3dDevice.Get(),
			1, (UINT)mAllRitems.size(), (UINT)mMaterials.size()));
	}
}

void LitColumnsApp::BuildMaterials()
{
	int i = 0;
	auto bricks0 = std::make_unique<Material>();
	bricks0->Name = "bricks0";
	bricks0->MatCBIndex = i;
	bricks0->DiffuseSrvHeapIndex = i++;
	bricks0->DiffuseAlbedo = XMFLOAT4(Colors::ForestGreen);
	bricks0->FresnelR0 = XMFLOAT3(0.02f, 0.02f, 0.02f);
	bricks0->Roughness = 0.1f;

	auto stone0 = std::make_unique<Material>();
	stone0->Name = "stone0";
	stone0->MatCBIndex = i;
	stone0->DiffuseSrvHeapIndex = i++;
	stone0->DiffuseAlbedo = XMFLOAT4(Colors::LightSteelBlue);
	stone0->FresnelR0 = XMFLOAT3(0.05f, 0.05f, 0.05f);
	stone0->Roughness = 0.3f;

	auto tile0 = std::make_unique<Material>();
	tile0->Name = "tile0";
	tile0->MatCBIndex = i;
	tile0->DiffuseSrvHeapIndex = i++;
	tile0->DiffuseAlbedo = XMFLOAT4(Colors::LightGray);
	tile0->FresnelR0 = XMFLOAT3(0.02f, 0.02f, 0.02f);
	tile0->Roughness = 0.2f;


	auto water0 = std::make_unique<Material>();
	water0->Name = "water0";
	water0->MatCBIndex = i;
	water0->DiffuseSrvHeapIndex = i++;
	water0->DiffuseAlbedo = XMFLOAT4(Colors::Blue);
	water0->FresnelR0 = XMFLOAT3(0.02f, 0.02f, 0.02f);
	water0->Roughness = 0.1f;

	auto grass0 = std::make_unique<Material>();
	grass0->Name = "grass0";
	grass0->MatCBIndex = i;
	grass0->DiffuseSrvHeapIndex = i++;
	grass0->DiffuseAlbedo = XMFLOAT4(Colors::LightGreen);
	grass0->FresnelR0 = XMFLOAT3(0.02f, 0.02f, 0.02f);
	grass0->Roughness = 0.1f;

	auto wood0 = std::make_unique<Material>();
	wood0->Name = "wood0";
	wood0->MatCBIndex = i;
	wood0->DiffuseSrvHeapIndex = i++;
	wood0->DiffuseAlbedo = XMFLOAT4(Colors::Beige);
	wood0->FresnelR0 = XMFLOAT3(0.02f, 0.02f, 0.02f);
	wood0->Roughness = 0.1f;

	auto grey0 = std::make_unique<Material>();
	grey0->Name = "grey0";
	grey0->MatCBIndex = i;
	grey0->DiffuseSrvHeapIndex = i++;
	grey0->DiffuseAlbedo = XMFLOAT4(Colors::LightGray);
	grey0->FresnelR0 = XMFLOAT3(0.02f, 0.02f, 0.02f);
	grey0->Roughness = 0.1f;

	auto black0 = std::make_unique<Material>();
	black0->Name = "black0";
	black0->MatCBIndex = i;
	black0->DiffuseSrvHeapIndex = i++;
	black0->DiffuseAlbedo = XMFLOAT4(Colors::Black);
	black0->FresnelR0 = XMFLOAT3(0.02f, 0.02f, 0.02f);
	black0->Roughness = 0.1f;

	auto darkGrey0 = std::make_unique<Material>();
	darkGrey0->Name = "darkGrey0";
	darkGrey0->MatCBIndex = i;
	darkGrey0->DiffuseSrvHeapIndex = i++;
	darkGrey0->DiffuseAlbedo = XMFLOAT4(Colors::Gray);
	darkGrey0->FresnelR0 = XMFLOAT3(0.02f, 0.02f, 0.02f);
	darkGrey0->Roughness = 0.1f;

	auto darkerGrey0 = std::make_unique<Material>();
	darkerGrey0->Name = "darkerGrey0";
	darkerGrey0->MatCBIndex = i;
	darkerGrey0->DiffuseSrvHeapIndex = i++;
	darkerGrey0->DiffuseAlbedo = XMFLOAT4(Colors::DarkSlateGray);
	darkerGrey0->FresnelR0 = XMFLOAT3(0.02f, 0.02f, 0.02f);
	darkerGrey0->Roughness = 0.1f;

	auto bush0 = std::make_unique<Material>();
	bush0->Name = "bush0";
	bush0->MatCBIndex = i;
	bush0->DiffuseSrvHeapIndex = i++;
	bush0->DiffuseAlbedo = XMFLOAT4(Colors::DarkGreen);
	bush0->FresnelR0 = XMFLOAT3(0.02f, 0.02f, 0.02f);
	bush0->Roughness = 0.1f;

	auto gold0 = std::make_unique<Material>();
	gold0->Name = "gold0";
	gold0->MatCBIndex = i;
	gold0->DiffuseSrvHeapIndex = i++;
	gold0->DiffuseAlbedo = XMFLOAT4(Colors::Gold);
	gold0->FresnelR0 = XMFLOAT3(0.02f, 0.02f, 0.02f);
	gold0->Roughness = 0.1f;

	auto red0 = std::make_unique<Material>();
	red0->Name = "red0";
	red0->MatCBIndex = i;
	red0->DiffuseSrvHeapIndex = i++;
	red0->DiffuseAlbedo = XMFLOAT4(Colors::Red);
	red0->FresnelR0 = XMFLOAT3(0.02f, 0.02f, 0.02f);
	red0->Roughness = 0.1f;

	mMaterials["bricks0"] = std::move(bricks0);
	mMaterials["stone0"] = std::move(stone0);
	mMaterials["tile0"] = std::move(tile0);
	mMaterials["water0"] = std::move(water0);
	mMaterials["grass0"] = std::move(grass0);
	mMaterials["wood0"] = std::move(wood0);
	mMaterials["grey0"] = std::move(grey0);
	mMaterials["black0"] = std::move(black0);
	mMaterials["darkGrey0"] = std::move(darkGrey0);
	mMaterials["darkerGrey0"] = std::move(darkerGrey0);
	mMaterials["bush0"] = std::move(bush0);
	mMaterials["gold0"] = std::move(gold0);
	mMaterials["red0"] = std::move(red0);

}

void LitColumnsApp::BuildRenderItem(int index,  std::string itemName, std::string material, XMFLOAT3 scaling, XMFLOAT3 rotation, XMFLOAT3 translation) {

	auto genItem = std::make_unique<RenderItem>();

	XMMATRIX transResult = XMMatrixScaling(scaling.x, scaling.y, scaling.z);
	transResult *= XMMatrixRotationRollPitchYaw(XMConvertToRadians(rotation.x), XMConvertToRadians(rotation.y), XMConvertToRadians(rotation.z));
	transResult *= XMMatrixTranslation(translation.x, translation.y, translation.z);
	XMStoreFloat4x4(&genItem->World, transResult);
	genItem->ObjCBIndex = index;
	genItem->Mat = mMaterials[material].get();
	genItem->Geo = mGeometries["shapeGeo"].get();
	genItem->PrimitiveType = D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST;
	genItem->IndexCount = genItem->Geo->DrawArgs[itemName].IndexCount;
	genItem->StartIndexLocation = genItem->Geo->DrawArgs[itemName].StartIndexLocation;
	genItem->BaseVertexLocation = genItem->Geo->DrawArgs[itemName].BaseVertexLocation;
	mAllRitems.push_back(std::move(genItem));
}

void LitColumnsApp::BuildRenderItems()
{
	int i = 0;
	//// Water
	BuildRenderItem(i++, "cube", "water0", XMFLOAT3(75, 9.26f, 75), XMFLOAT3(0, 0, 0), XMFLOAT3(0, -8.4f, 0));
	//// Grass
	BuildRenderItem(i++, "cube", "grass0", XMFLOAT3(50, 9.26f, 50), XMFLOAT3(0, 0, 0), XMFLOAT3(0, -4.53f, 0));
	BuildRenderItem(i++, "cube", "grass0", XMFLOAT3(150, 9.26f, 50), XMFLOAT3(0, 0, 0), XMFLOAT3(0, -4.53f, 60.3f));
	BuildRenderItem(i++, "cube", "grass0", XMFLOAT3(150, 9.26, 50), XMFLOAT3(0, 0, 0), XMFLOAT3(0, -4.53f, -59.2f));
	BuildRenderItem(i++, "cube", "grass0", XMFLOAT3(38.52, 9.26, 100), XMFLOAT3(0, 0, 0), XMFLOAT3(55.74, -4.53, 0));
	BuildRenderItem(i++, "cube", "grass0", XMFLOAT3(38.52, 9.26, 100), XMFLOAT3(0, 0, 0), XMFLOAT3(-55.56, -4.53, 0));
	//// Bridge
	BuildRenderItem(i++, "cube", "wood0", XMFLOAT3(19.05, 1, 8.33), XMFLOAT3(0, 0, 0), XMFLOAT3(30.24, 0.39, 0));
	//// Walls
	BuildRenderItem(i++, "cube", "grey0", XMFLOAT3(27.86, 1.76, 44.744), XMFLOAT3(0, 180, 90), XMFLOAT3(-22.27, 2, 0.144));
	BuildRenderItem(i++, "cube", "grey0", XMFLOAT3(27.86, 1.76, 44.744), XMFLOAT3(0, 0, 90), XMFLOAT3(20.74, 2, 0.144));
	BuildRenderItem(i++, "cube", "grey0", XMFLOAT3(27.86, 1.76, 44.744), XMFLOAT3(0, 90, 90), XMFLOAT3(-0.72, 2, -21.68));
	BuildRenderItem(i++, "cube", "grey0", XMFLOAT3(27.86, 1.76, 44.744), XMFLOAT3(0, 270, 90), XMFLOAT3(-0.72, 2, 21.49));
	//// Central Cube
	BuildRenderItem(i++, "cube", "grey0", XMFLOAT3(23, 15, 23), XMFLOAT3(0, 180, 90), XMFLOAT3(0, 4.4, 0));
	//// Portal
	BuildRenderItem(i++, "cube", "black0", XMFLOAT3(19.05, 1.95, 8.32), XMFLOAT3(0, 0, 0), XMFLOAT3(20.76, -2.4, 0));
	BuildRenderItem(i++, "myCylinder", "black0", XMFLOAT3(8.38, 0.96, 8.38), XMFLOAT3(0, 0, 0), XMFLOAT3(20.76, 7.02, 0.011));
	//// Towers
	BuildRenderItem(i++, "myCylinder", "darkGrey0", XMFLOAT3(4, 40, 4), XMFLOAT3(0, 0, 0), XMFLOAT3(20.8, 2.73, -21.5));
	BuildRenderItem(i++, "myCylinder", "darkGrey0", XMFLOAT3(4, 40, 4), XMFLOAT3(0, 0, 0), XMFLOAT3(-22.21, 2.73, -21.5));
	BuildRenderItem(i++, "myCylinder", "darkGrey0", XMFLOAT3(4, 40, 4), XMFLOAT3(0, 0, 0), XMFLOAT3(20.8, 2.73, 21.67));
	BuildRenderItem(i++, "myCylinder", "darkGrey0", XMFLOAT3(4, 40, 4), XMFLOAT3(0, 0, 0), XMFLOAT3(-20.19, 2.73, 21.67));
	BuildRenderItem(i++, "myCylinder", "darkerGrey0", XMFLOAT3(5, 0.3, 5), XMFLOAT3(0, 0, 0), XMFLOAT3(20.8, 22.83, -21.5));
	BuildRenderItem(i++, "myCylinder", "darkerGrey0", XMFLOAT3(5, 0.3, 5), XMFLOAT3(0, 0, 0), XMFLOAT3(-22.21, 22.83, -21.5));
	BuildRenderItem(i++, "myCylinder", "darkerGrey0", XMFLOAT3(5, 0.3, 5), XMFLOAT3(0, 0, 0), XMFLOAT3(20.8, 22.83, 21.67));
	BuildRenderItem(i++, "myCylinder", "darkerGrey0", XMFLOAT3(5, 0.3, 5), XMFLOAT3(0, 0, 0), XMFLOAT3(-20.19, 22.83, 21.67));
	//// Trees
	BuildRenderItem(i++, "myCylinder", "wood0", XMFLOAT3(1, 10, 1), XMFLOAT3(0, 0, 0), XMFLOAT3(45, 2, -20));
	BuildRenderItem(i++, "myCylinder", "wood0", XMFLOAT3(1, 10, 1), XMFLOAT3(0, 0, 0), XMFLOAT3(45, 2, -10));
	BuildRenderItem(i++, "myCylinder", "wood0", XMFLOAT3(1, 10, 1), XMFLOAT3(0, 0, 0), XMFLOAT3(45, 2, 0));
	BuildRenderItem(i++, "myCylinder", "wood0", XMFLOAT3(1, 10, 1), XMFLOAT3(0, 0, 0), XMFLOAT3(45, 2, 10));
	BuildRenderItem(i++, "myCylinder", "wood0", XMFLOAT3(1, 10, 1), XMFLOAT3(0, 0, 0), XMFLOAT3(45, 2, 20));
	//// Treetops
	BuildRenderItem(i++, "cone", "bush0", XMFLOAT3(5, 10, 5), XMFLOAT3(0, 0, 0), XMFLOAT3(45, 10, -20));
	BuildRenderItem(i++, "cone", "bush0", XMFLOAT3(5, 10, 5), XMFLOAT3(0, 0, 0), XMFLOAT3(45, 10, -10));
	BuildRenderItem(i++, "cone", "bush0", XMFLOAT3(5, 10, 5), XMFLOAT3(0, 0, 0), XMFLOAT3(45, 10, 0));
	BuildRenderItem(i++, "cone", "bush0", XMFLOAT3(5, 10, 5), XMFLOAT3(0, 0, 0), XMFLOAT3(45, 10, 10));
	BuildRenderItem(i++, "cone", "bush0", XMFLOAT3(5, 10, 5), XMFLOAT3(0, 0, 0), XMFLOAT3(45, 10, 20));
	//// Strange sphere on top of central building
	BuildRenderItem(i++, "mySphere", "darkGrey0", XMFLOAT3(15, 15, 15), XMFLOAT3(0, 0, 0), XMFLOAT3(0, 11.91, 0));
	//// Strange indicator
	BuildRenderItem(i++, "simsInd", "gold0", XMFLOAT3(1,1,1), XMFLOAT3(90, 0, 0), XMFLOAT3(21.55, 5.9, -11));
	BuildRenderItem(i++, "simsInd", "gold0", XMFLOAT3(1,1,1), XMFLOAT3(90, 0, 0), XMFLOAT3(21.55, 5.9, 11));
	//// StarPortal
	BuildRenderItem(i++, "star", "gold0", XMFLOAT3(15, 15, 1), XMFLOAT3(0, 90, 0), XMFLOAT3(21.55, 5.9, 0));
	//// Pyramids for good luck
	BuildRenderItem(i++, "pyramid", "gold0", XMFLOAT3(3, 5, 3), XMFLOAT3(0, 0, 0), XMFLOAT3(20.8, 24.83, -21.5));
	BuildRenderItem(i++, "pyramid", "gold0", XMFLOAT3(3, 5, 3), XMFLOAT3(0, 0, 0), XMFLOAT3(-22.21, 24.83, -21.5));
	BuildRenderItem(i++, "pyramid", "gold0", XMFLOAT3(3, 5, 3), XMFLOAT3(0, 0, 0), XMFLOAT3(20.8, 24.83, 21.67));
	BuildRenderItem(i++, "pyramid", "gold0", XMFLOAT3(3, 5, 3), XMFLOAT3(0, 0, 0), XMFLOAT3(-20.19, 24.83, 21.67));
	//// The spaceship
	BuildRenderItem(i++, "hexagon", "gold0", XMFLOAT3(50, 5, 50), XMFLOAT3(0, 0, 0), XMFLOAT3(0, 40, 0));
	//// Pathways
	BuildRenderItem(i++, "wedge", "red0", XMFLOAT3(2, 2, 2), XMFLOAT3(0, 90, 0), XMFLOAT3(39.35, 0.26, 5.51));
	BuildRenderItem(i++, "wedge", "red0", XMFLOAT3(2, 2, 2), XMFLOAT3(0, 90, 0), XMFLOAT3(41.35, 0.26, 5.51));
	BuildRenderItem(i++, "wedge", "red0", XMFLOAT3(2, 2, 2), XMFLOAT3(0, 90, 0), XMFLOAT3(43.35, 0.26, 5.51));
	BuildRenderItem(i++, "wedge", "red0", XMFLOAT3(2, 2, 2), XMFLOAT3(0, 90, 0), XMFLOAT3(45.35, 0.26, 5.51));
	BuildRenderItem(i++, "wedge", "red0", XMFLOAT3(2, 2, 2), XMFLOAT3(0, 90, 0), XMFLOAT3(47.35, 0.26, 5.51));
	BuildRenderItem(i++, "wedge", "red0", XMFLOAT3(2, 2, 2), XMFLOAT3(0, 90, 0), XMFLOAT3(39.35, 0.26, -2.58));
	BuildRenderItem(i++, "wedge", "red0", XMFLOAT3(2, 2, 2), XMFLOAT3(0, 90, 0), XMFLOAT3(41.35, 0.26, -2.58));
	BuildRenderItem(i++, "wedge", "red0", XMFLOAT3(2, 2, 2), XMFLOAT3(0, 90, 0), XMFLOAT3(43.35, 0.26, -2.58));
	BuildRenderItem(i++, "wedge", "red0", XMFLOAT3(2, 2, 2), XMFLOAT3(0, 90, 0), XMFLOAT3(45.35, 0.26, -2.58));
	BuildRenderItem(i++, "wedge", "red0", XMFLOAT3(2, 2, 2), XMFLOAT3(0, 90, 0), XMFLOAT3(47.35, 0.26, -2.58));
	//// Why a prism here?
	BuildRenderItem(i++, "triPrism", "red0", XMFLOAT3(5, 5, 5), XMFLOAT3(0, 0, 0), XMFLOAT3(0, 20, 0));

	// All the render items are opaque.
	for (auto& e : mAllRitems)
		mOpaqueRitems.push_back(e.get());
}

void LitColumnsApp::DrawRenderItems(ID3D12GraphicsCommandList* cmdList, const std::vector<RenderItem*>& ritems)
{
	UINT objCBByteSize = d3dUtil::CalcConstantBufferByteSize(sizeof(ObjectConstants));
	UINT matCBByteSize = d3dUtil::CalcConstantBufferByteSize(sizeof(MaterialConstants));

	auto objectCB = mCurrFrameResource->ObjectCB->Resource();
	auto matCB = mCurrFrameResource->MaterialCB->Resource();

	// For each render item...
	for (size_t i = 0; i < ritems.size(); ++i)
	{
		auto ri = ritems[i];

		cmdList->IASetVertexBuffers(0, 1, &ri->Geo->VertexBufferView());
		cmdList->IASetIndexBuffer(&ri->Geo->IndexBufferView());
		cmdList->IASetPrimitiveTopology(ri->PrimitiveType);

		D3D12_GPU_VIRTUAL_ADDRESS objCBAddress = objectCB->GetGPUVirtualAddress() + ri->ObjCBIndex*objCBByteSize;
		D3D12_GPU_VIRTUAL_ADDRESS matCBAddress = matCB->GetGPUVirtualAddress() + ri->Mat->MatCBIndex*matCBByteSize;

		cmdList->SetGraphicsRootConstantBufferView(0, objCBAddress);
		cmdList->SetGraphicsRootConstantBufferView(1, matCBAddress);

		cmdList->DrawIndexedInstanced(ri->IndexCount, 1, ri->StartIndexLocation, ri->BaseVertexLocation, 0);
	}
}
