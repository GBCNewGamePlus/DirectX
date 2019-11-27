#include "../../Common/d3dApp.h"
#include "../../Common/MathHelper.h"
#include "../../Common/UploadBuffer.h"
#include "GeometryGenerator.h"
#include "FrameResource.h"
#include "Waves.h"
#include "Camera.h"
#include <tchar.h>

using Microsoft::WRL::ComPtr;
using namespace DirectX;
using namespace DirectX::PackedVector;
using namespace std;

#pragma comment(lib, "d3dcompiler.lib")
#pragma comment(lib, "D3D12.lib")

const int gNumFrameResources = 3;

struct BlockInfo
{
	XMFLOAT3 position;
	XMFLOAT3 scale;

	BlockInfo()
	{
		position = XMFLOAT3(0, 0, 0);
		scale = XMFLOAT3(0, 0, 0);
	}

	BlockInfo(XMFLOAT3 pos, XMFLOAT3 sca)
	{
		position = pos;
		scale = sca;
	}
};

BlockInfo blocks[1];
int blocksSize = 1;
int collisionOffset = 2;
bool blockX, blockZ;

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
	D3D12_PRIMITIVE_TOPOLOGY PrimitiveType = D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST;

	// DrawIndexedInstanced parameters.
	UINT IndexCount = 0;
	UINT StartIndexLocation = 0;
	int BaseVertexLocation = 0;
};

enum class RenderLayer : int
{
	Opaque = 0,
	Transparent,
	AlphaTested,
	AlphaTestedTreeSprites,
	Count
};

class Castle : public D3DApp
{
public:
	Castle(HINSTANCE hInstance);
	Castle(const Castle& rhs) = delete;
	Castle& operator=(const Castle& rhs) = delete;
	~Castle();

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
	void UpdateWaves(const GameTimer& gt);

	void LoadTextures();
	void BuildRootSignature();
	void BuildDescriptorHeaps();
	void BuildShadersAndInputLayouts();
	void BuildLandGeometry();
	void BuildWavesGeometry();
	void BuildBoxGeometry();
	void BuildTreeSpritesGeometry();
	void BuildShapeGeometry();
	void BuildPSOs();
	void BuildFrameResources();
	void BuildMaterials();
	void BuildRenderGeoItems();
	void DrawRenderItems(ID3D12GraphicsCommandList* cmdList, const std::vector<RenderItem*>& ritems);
	void CollisionDetection();

	std::array<const CD3DX12_STATIC_SAMPLER_DESC, 6> GetStaticSamplers();

	float GetHillsHeight(float x, float z)const;
	XMFLOAT3 GetHillsNormal(float x, float z)const;

	void BuildRenderGeoItem(int index, std::string itemName, std::string material, XMFLOAT3 scaling, XMFLOAT3 rotation, XMFLOAT3 translation);

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
	std::unordered_map<std::string, ComPtr<ID3D12PipelineState>> mPSOs;

	std::vector<D3D12_INPUT_ELEMENT_DESC> mStdInputLayout;
	std::vector<D3D12_INPUT_ELEMENT_DESC> mTreeSpriteInputLayout;

	RenderItem* mWavesRitem = nullptr;

	// List of all the render items.
	std::vector<std::unique_ptr<RenderItem>> mAllRitems;

	// Render items divided by PSO.
	std::vector<RenderItem*> mRitemLayer[(int)RenderLayer::Count];

	std::unique_ptr<Waves> mWaves;

	PassConstants mMainPassCB;
	Camera mCamera;

	POINT mLastMousePos;

	float mBasePathway = 42;
	float mPathwayLeft = 4.58f;
	float mPathwayRight = -4.58f;

};

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE prevInstance, PSTR cmdLine, int showCmd)
{
	// Enable run-time memory check for debug builds.
#if defined(DEBUG) | defined(_DEBUG)
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif

	try
	{
		Castle theApp(hInstance);
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

Castle::Castle(HINSTANCE hInstance)
	: D3DApp(hInstance)
{
}

Castle::~Castle()
{
	if (md3dDevice != nullptr)
		FlushCommandQueue();
}

bool Castle::Initialize()
{
	if (!D3DApp::Initialize())
		return false;

	// Reset the command list to prep for initialization commands.
	ThrowIfFailed(mCommandList->Reset(mDirectCmdListAlloc.Get(), nullptr));

	// Get the increment size of a descriptor in this heap type.  This is hardware specific, 
	// so we have to query this information.
	mCbvSrvDescriptorSize = md3dDevice->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV);

	mWaves = std::make_unique<Waves>(155, 155, 1.0f, 0.03f, 4.0f, 0.2f);

	mCamera.SetPosition(67.0f, mCamera.YPOSITION, -1.2f);

	LoadTextures();
	BuildRootSignature();
	BuildDescriptorHeaps();
	BuildShadersAndInputLayouts();
	BuildLandGeometry();
	BuildWavesGeometry();
	BuildBoxGeometry();
	BuildTreeSpritesGeometry();
	BuildShapeGeometry();
	BuildMaterials();
	BuildRenderGeoItems();
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

void Castle::OnResize()
{
	D3DApp::OnResize();
	// The window resized, so update the aspect ratio and recompute the projection matrix.
	mCamera.SetLens(0.25f*MathHelper::Pi, AspectRatio(), 1.0f, 1000.0f);
}

void Castle::Update(const GameTimer& gt)
{
	CollisionDetection();
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
	UpdateWaves(gt);
}

void Castle::CollisionDetection()
{
	blockX = false;
	blockZ = false;
	float cameraX = mCamera.GetPosition3f().x;
	float cameraZ = mCamera.GetPosition3f().z;

	for (int i = 0; i < blocksSize; i++)
	{
		float cubeX1 = blocks[i].position.x - blocks[i].scale.x / 2;
		float cubeX2 = blocks[i].position.x + blocks[i].scale.x / 2;
		float cubeZ1 = blocks[i].position.z - blocks[i].scale.z / 2;
		float cubeZ2 = blocks[i].position.z + blocks[i].scale.z / 2;

		char buf[100];
		sprintf_s(buf, "X:%f - %f\nZ:%f - %f\nCamX:%f CamZ:%f\n\n", cubeX1, cubeX2, cubeZ1, cubeZ2, cameraX, cameraZ);
		OutputDebugStringA(buf);
		if (cameraX >= cubeX1 && cameraX <= cubeX2)
		{
			if (cameraZ >= cubeZ1 && cameraZ <= cubeZ2)
			{
				mCamera.SetPosition(0, 0, 0);
			}
		}
	}
}

void Castle::Draw(const GameTimer& gt)
{
	auto cmdListAlloc = mCurrFrameResource->CmdListAlloc;

	// Reuse the memory associated with command recording.
	// We can only reset when the associated command lists have finished execution on the GPU.
	ThrowIfFailed(cmdListAlloc->Reset());

	// A command list can be reset after it has been added to the command queue via ExecuteCommandList.
	// Reusing the command list reuses memory.
	ThrowIfFailed(mCommandList->Reset(cmdListAlloc.Get(), mPSOs["opaque"].Get()));

	mCommandList->RSSetViewports(1, &mScreenViewport);
	mCommandList->RSSetScissorRects(1, &mScissorRect);

	// Indicate a state transition on the resource usage.
	mCommandList->ResourceBarrier(1, &CD3DX12_RESOURCE_BARRIER::Transition(CurrentBackBuffer(),
		D3D12_RESOURCE_STATE_PRESENT, D3D12_RESOURCE_STATE_RENDER_TARGET));

	// Clear the back buffer and depth buffer.
	mCommandList->ClearRenderTargetView(CurrentBackBufferView(), (float*)&mMainPassCB.FogColor, 0, nullptr);
	mCommandList->ClearDepthStencilView(DepthStencilView(), D3D12_CLEAR_FLAG_DEPTH | D3D12_CLEAR_FLAG_STENCIL, 1.0f, 0, 0, nullptr);

	// Specify the buffers we are going to render to.
	mCommandList->OMSetRenderTargets(1, &CurrentBackBufferView(), true, &DepthStencilView());

	ID3D12DescriptorHeap* descriptorHeaps[] = { mSrvDescriptorHeap.Get() };
	mCommandList->SetDescriptorHeaps(_countof(descriptorHeaps), descriptorHeaps);

	mCommandList->SetGraphicsRootSignature(mRootSignature.Get());

	auto passCB = mCurrFrameResource->PassCB->Resource();
	mCommandList->SetGraphicsRootConstantBufferView(2, passCB->GetGPUVirtualAddress());

	DrawRenderItems(mCommandList.Get(), mRitemLayer[(int)RenderLayer::Opaque]);

	mCommandList->SetPipelineState(mPSOs["alphaTested"].Get());
	DrawRenderItems(mCommandList.Get(), mRitemLayer[(int)RenderLayer::AlphaTested]);

	mCommandList->SetPipelineState(mPSOs["treeSprites"].Get());
	DrawRenderItems(mCommandList.Get(), mRitemLayer[(int)RenderLayer::AlphaTestedTreeSprites]);

	mCommandList->SetPipelineState(mPSOs["transparent"].Get());
	DrawRenderItems(mCommandList.Get(), mRitemLayer[(int)RenderLayer::Transparent]);

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

void Castle::OnMouseDown(WPARAM btnState, int x, int y)
{
	mLastMousePos.x = x;
	mLastMousePos.y = y;

	SetCapture(mhMainWnd);
}

void Castle::OnMouseUp(WPARAM btnState, int x, int y)
{
	ReleaseCapture();
}

void Castle::OnMouseMove(WPARAM btnState, int x, int y)
{
	//if left mouse button is down and moving
	if ((btnState & MK_LBUTTON) != 0)
	{
		// Make each pixel correspond to a quarter of a degree.
		float dx = XMConvertToRadians(0.25f*static_cast<float>(x - mLastMousePos.x));
		float dy = XMConvertToRadians(0.25f*static_cast<float>(y - mLastMousePos.y));

		mCamera.Pitch(dy);
		mCamera.RotateY(dx);
	}

	mLastMousePos.x = x;
	mLastMousePos.y = y;
}

void Castle::OnKeyboardInput(const GameTimer& gt)
{
	const float dt = gt.DeltaTime();

	//GetAsyncKeyState returns a short (2 bytes)
	if (GetAsyncKeyState('W') & 0x8000) //most significant bit (MSB) is 1 when key is pressed (1000 000 000 000)
		mCamera.Walk(10.0f*dt);

	if (GetAsyncKeyState('S') & 0x8000)
		mCamera.Walk(-10.0f*dt);

	if (GetAsyncKeyState('A') & 0x8000)
		mCamera.Strafe(-10.0f*dt);

	if (GetAsyncKeyState('D') & 0x8000)
		mCamera.Strafe(10.0f*dt);

	mCamera.UpdateViewMatrix();
}

void Castle::UpdateCamera(const GameTimer& gt)
{
	// Convert Spherical to Cartesian coordinates.
	//mEyePos.x = 67.0f;// mRadius * sinf(mPhi)*cosf(mTheta); //67.0f
	//mEyePos.z = -1.2f;// mRadius * sinf(mPhi)*sinf(mTheta); //-1.2f
	//mEyePos.y = 5.5f;// 

	//// Build the view matrix.
	//XMVECTOR pos = XMVectorSet(mEyePos.x, mEyePos.y, mEyePos.z, 1.0f);
	//XMVECTOR target = XMVectorZero();
	//XMVECTOR up = XMVectorSet(0.0f, 1.0f, 0.0f, 0.0f);

	//XMMATRIX view = XMMatrixLookAtLH(pos, target, up);
	//XMStoreFloat4x4(&mView, view);
}

void Castle::AnimateMaterials(const GameTimer& gt)
{
	// Scroll the water material texture coordinates.
	auto waterMat = mMaterials["water"].get();

	float& tu = waterMat->MatTransform(3, 0);
	float& tv = waterMat->MatTransform(3, 1);

	tu += 0.1f * gt.DeltaTime();
	tv += 0.02f * gt.DeltaTime();

	if (tu >= 1.0f)
		tu -= 1.0f;

	if (tv >= 1.0f)
		tv -= 1.0f;

	waterMat->MatTransform(3, 0) = tu;
	waterMat->MatTransform(3, 1) = tv;

	// Material has changed, so need to update cbuffer.
	waterMat->NumFramesDirty = gNumFrameResources;
}

void Castle::UpdateObjectCBs(const GameTimer& gt)
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

void Castle::UpdateMaterialCBs(const GameTimer& gt)
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

void Castle::UpdateMainPassCB(const GameTimer& gt)
{

	XMMATRIX view = mCamera.GetView();
	XMMATRIX proj = mCamera.GetProj();

	//XMMATRIX view = XMLoadFloat4x4(&mView);
	//XMMATRIX proj = XMLoadFloat4x4(&mProj);

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

	mMainPassCB.EyePosW = mCamera.GetPosition3f();
	mMainPassCB.RenderTargetSize = XMFLOAT2((float)mClientWidth, (float)mClientHeight);
	mMainPassCB.InvRenderTargetSize = XMFLOAT2(1.0f / mClientWidth, 1.0f / mClientHeight);

	mMainPassCB.NearZ = 1.0f;
	mMainPassCB.FarZ = 1000.0f;
	mMainPassCB.TotalTime = gt.TotalTime();
	mMainPassCB.DeltaTime = gt.DeltaTime();
	mMainPassCB.AmbientLight = {1.0f, 1.0f, 1.0f, 1.0f };

	// SIMS IND 1 and 2
	mMainPassCB.Lights[0].Direction = { 1, 0, 0 };
	mMainPassCB.Lights[0].Strength = {0.8f, 0.1f, 0.1f};
	mMainPassCB.Lights[0].Position = {21.55f, 5.9f, -11};
	mMainPassCB.Lights[0].FalloffStart = 10;
	mMainPassCB.Lights[0].FalloffEnd = 12;

	mMainPassCB.Lights[1].Direction = { 1, 0, 0 };
	mMainPassCB.Lights[1].Strength = { 0.8f, 0.1f, 0.1f };
	mMainPassCB.Lights[1].Position = { 21.55f, 5.9f, 11 };
	mMainPassCB.Lights[1].FalloffStart = 10;
	mMainPassCB.Lights[1].FalloffEnd = 12;

	mMainPassCB.Lights[2].Direction = { 0, -1, 0 };
	mMainPassCB.Lights[2].Strength = { 1, 1 , 1 };
	mMainPassCB.Lights[2].Position = { 0, 50, 0 };
	mMainPassCB.Lights[2].FalloffStart = 15;
	mMainPassCB.Lights[2].FalloffEnd = 40;

	mMainPassCB.Lights[3].Direction = { 0, -1, 0 };
	mMainPassCB.Lights[3].Strength = { 0.8f, 0.1f, 0.1f };
	mMainPassCB.Lights[3].Position = { 20.8f, 26.83f, -21.5f };
	mMainPassCB.Lights[3].FalloffStart = 10;
	mMainPassCB.Lights[3].FalloffEnd = 12;

	mMainPassCB.Lights[4].Direction = { 0, -1, 0 };
	mMainPassCB.Lights[4].Strength = { 0.8f, 0.1f, 0.1f };
	mMainPassCB.Lights[4].Position = { -22.21f, 26.83f, -21.5f };
	mMainPassCB.Lights[4].FalloffStart = 10;
	mMainPassCB.Lights[4].FalloffEnd = 12;

	mMainPassCB.Lights[5].Direction = { 0, -1, 0 };
	mMainPassCB.Lights[5].Strength = { 0.8f, 0.1f, 0.1f };
	mMainPassCB.Lights[5].Position = { 20.8f, 26.83f, 21.67f };
	mMainPassCB.Lights[5].FalloffStart = 10;
	mMainPassCB.Lights[5].FalloffEnd = 12;

	mMainPassCB.Lights[6].Direction = { 0, -1, 0 };
	mMainPassCB.Lights[6].Strength = { 0.8f, 0.1f, 0.1f };
	mMainPassCB.Lights[6].Position = { -20.19f, 26.83f, 21.67f };
	mMainPassCB.Lights[6].FalloffStart = 10;
	mMainPassCB.Lights[6].FalloffEnd = 12;

	// Beacons
	mMainPassCB.Lights[8].Direction = { 0, -1, 0 };
	mMainPassCB.Lights[8].Strength = { 0.8f, 0.8f, 0.8f };
	mMainPassCB.Lights[8].Position = { mBasePathway + 10, 12.26f, mPathwayRight - 1 };
	mMainPassCB.Lights[8].FalloffStart = 10;
	mMainPassCB.Lights[8].FalloffEnd = 20;

	mMainPassCB.Lights[9].Direction = { 0, -1, 0 };
	mMainPassCB.Lights[9].Strength = { 0.8f, 0.8f, 0.8f };
	mMainPassCB.Lights[9].Position = { mBasePathway + 10, 12.26f, mPathwayLeft + 1 };
	mMainPassCB.Lights[9].FalloffStart = 10;
	mMainPassCB.Lights[9].FalloffEnd = 20;


	mMainPassCB.Lights[10].Direction = { -1, -1, 0 };
	mMainPassCB.Lights[10].Strength = { 0.8f, 0.8f, 0.8f };
	mMainPassCB.Lights[10].Position = { 23.55f, 15.9f, 0 };
	mMainPassCB.Lights[10].FalloffStart = 10;
	mMainPassCB.Lights[10].FalloffEnd = 40;

	auto currPassCB = mCurrFrameResource->PassCB.get();
	currPassCB->CopyData(0, mMainPassCB);
}

void Castle::UpdateWaves(const GameTimer& gt)
{
	// Every quarter second, generate a random wave.
	static float t_base = 0.0f;
	if ((mTimer.TotalTime() - t_base) >= 0.25f)
	{
		t_base += 0.25f;

		int i = MathHelper::Rand(4, mWaves->RowCount() - 5);
		int j = MathHelper::Rand(4, mWaves->ColumnCount() - 5);

		float r = MathHelper::RandF(0.1f, 0.3f);

		mWaves->Disturb(i, j, r);
	}

	// Update the wave simulation.
	mWaves->Update(gt.DeltaTime());

	// Update the wave vertex buffer with the new solution.
	auto currWavesVB = mCurrFrameResource->WavesVB.get();
	for (int i = 0; i < mWaves->VertexCount(); ++i)
	{
		Vertex v;

		v.Pos = mWaves->Position(i);
		v.Normal = mWaves->Normal(i);

		// Derive tex-coords from position by 
		// mapping [-w/2,w/2] --> [0,1]
		v.TexC.x = 0.5f + v.Pos.x / mWaves->Width();
		v.TexC.y = - v.Pos.z / mWaves->Depth();

		currWavesVB->CopyData(i, v);
	}

	// Set the dynamic VB of the wave renderitem to the current frame VB.
	mWavesRitem->Geo->VertexBufferGPU = currWavesVB->Resource();
}

void Castle::LoadTextures()
{
	auto grassTex = std::make_unique<Texture>();
	grassTex->Name = "grassTex";
	grassTex->Filename = L"../../Textures/grass.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(md3dDevice.Get(),
		mCommandList.Get(), grassTex->Filename.c_str(),
		grassTex->Resource, grassTex->UploadHeap));

	auto waterTex = std::make_unique<Texture>();
	waterTex->Name = "waterTex";
	waterTex->Filename = L"../../Textures/water1.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(md3dDevice.Get(),
		mCommandList.Get(), waterTex->Filename.c_str(),
		waterTex->Resource, waterTex->UploadHeap));

	auto fenceTex = std::make_unique<Texture>();
	fenceTex->Name = "fenceTex";
	fenceTex->Filename = L"../../Textures/WireFence.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(md3dDevice.Get(),
		mCommandList.Get(), fenceTex->Filename.c_str(),
		fenceTex->Resource, fenceTex->UploadHeap));

	auto treeArrayTex = std::make_unique<Texture>();
	treeArrayTex->Name = "treeArrayTex";
	treeArrayTex->Filename = L"../../Textures/treeArray2.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(md3dDevice.Get(),
		mCommandList.Get(), treeArrayTex->Filename.c_str(),
		treeArrayTex->Resource, treeArrayTex->UploadHeap));

	auto brickArrayTex = std::make_unique<Texture>();
	brickArrayTex->Name = "brickArrayTex";
	brickArrayTex->Filename = L"../../Textures/bricks2.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(md3dDevice.Get(),
		mCommandList.Get(), brickArrayTex->Filename.c_str(),
		brickArrayTex->Resource, brickArrayTex->UploadHeap));

	auto brickDarkTex = std::make_unique<Texture>();
	brickDarkTex->Name = "brickDarkTex";
	brickDarkTex->Filename = L"../../Textures/bricks4.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(md3dDevice.Get(),
		mCommandList.Get(), brickDarkTex->Filename.c_str(),
		brickDarkTex->Resource, brickDarkTex->UploadHeap));

	auto bricks3 = std::make_unique<Texture>();
	bricks3->Name = "bricks3";
	bricks3->Filename = L"../../Textures/bricks3.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(md3dDevice.Get(),
		mCommandList.Get(), bricks3->Filename.c_str(),
		bricks3->Resource, bricks3->UploadHeap));

	auto stone = std::make_unique<Texture>();
	stone->Name = "stone";
	stone->Filename = L"../../Textures/stone.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(md3dDevice.Get(),
		mCommandList.Get(), stone->Filename.c_str(),
		stone->Resource, stone->UploadHeap));

	auto black = std::make_unique<Texture>();
	black->Name = "black";
	black->Filename = L"../../Textures/black.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(md3dDevice.Get(),
		mCommandList.Get(), black->Filename.c_str(),
		black->Resource, black->UploadHeap));

	auto sims = std::make_unique<Texture>();
	sims->Name = "sims";
	sims->Filename = L"../../Textures/sims.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(md3dDevice.Get(),
		mCommandList.Get(), sims->Filename.c_str(),
		sims->Resource, sims->UploadHeap));

	auto bricksOnAcid = std::make_unique<Texture>();
	bricksOnAcid->Name = "bricksOnAcid";
	bricksOnAcid->Filename = L"../../Textures/bricksOnAcid.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(md3dDevice.Get(),
		mCommandList.Get(), bricksOnAcid->Filename.c_str(),
		bricksOnAcid->Resource, bricksOnAcid->UploadHeap));

	mTextures[grassTex->Name] = std::move(grassTex);
	mTextures[waterTex->Name] = std::move(waterTex);
	mTextures[fenceTex->Name] = std::move(fenceTex);
	mTextures[treeArrayTex->Name] = std::move(treeArrayTex);
	mTextures[brickArrayTex->Name] = std::move(brickArrayTex);
	mTextures[brickDarkTex->Name] = std::move(brickDarkTex);
	mTextures[bricks3->Name] = std::move(bricks3);
	mTextures[stone->Name] = std::move(stone);
	mTextures[black->Name] = std::move(black);
	mTextures[sims->Name] = std::move(sims);
	mTextures[bricksOnAcid->Name] = std::move(bricksOnAcid);
}

void Castle::BuildRootSignature()
{
	CD3DX12_DESCRIPTOR_RANGE texTable;
	texTable.Init(D3D12_DESCRIPTOR_RANGE_TYPE_SRV, 1, 0);

	// Root parameter can be a table, root descriptor or root constants.
	CD3DX12_ROOT_PARAMETER slotRootParameter[4];

	// Perfomance TIP: Order from most frequent to least frequent.
	slotRootParameter[0].InitAsDescriptorTable(1, &texTable, D3D12_SHADER_VISIBILITY_PIXEL);
	slotRootParameter[1].InitAsConstantBufferView(0);
	slotRootParameter[2].InitAsConstantBufferView(1);
	slotRootParameter[3].InitAsConstantBufferView(2);

	auto staticSamplers = GetStaticSamplers();

	// A root signature is an array of root parameters.
	CD3DX12_ROOT_SIGNATURE_DESC rootSigDesc(4, slotRootParameter,
		(UINT)staticSamplers.size(), staticSamplers.data(),
		D3D12_ROOT_SIGNATURE_FLAG_ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT);

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

void Castle::BuildDescriptorHeaps()
{
	//
	// Create the SRV heap.
	//
	D3D12_DESCRIPTOR_HEAP_DESC srvHeapDesc = {};
	srvHeapDesc.NumDescriptors = 11;
	srvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
	srvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;
	ThrowIfFailed(md3dDevice->CreateDescriptorHeap(&srvHeapDesc, IID_PPV_ARGS(&mSrvDescriptorHeap)));

	//
	// Fill out the heap with actual descriptors.
	//
	CD3DX12_CPU_DESCRIPTOR_HANDLE hDescriptor(mSrvDescriptorHeap->GetCPUDescriptorHandleForHeapStart());

	auto grassTex = mTextures["grassTex"]->Resource;
	auto waterTex = mTextures["waterTex"]->Resource;
	auto fenceTex = mTextures["fenceTex"]->Resource;
	auto treeArrayTex = mTextures["treeArrayTex"]->Resource;
	auto brickTex = mTextures["brickArrayTex"]->Resource;
	auto brickDarkTex = mTextures["brickDarkTex"]->Resource;
	auto bricks3 = mTextures["bricks3"]->Resource;
	auto stone = mTextures["stone"]->Resource;
	auto black = mTextures["black"]->Resource;
	auto sims = mTextures["sims"]->Resource;
	auto bricksOnAcid = mTextures["bricksOnAcid"]->Resource;

	D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
	srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
	srvDesc.Format = grassTex->GetDesc().Format;
	srvDesc.ViewDimension = D3D12_SRV_DIMENSION_TEXTURE2D;
	srvDesc.Texture2D.MostDetailedMip = 0;
	srvDesc.Texture2D.MipLevels = -1;
	md3dDevice->CreateShaderResourceView(grassTex.Get(), &srvDesc, hDescriptor);

	// next descriptor
	hDescriptor.Offset(1, mCbvSrvDescriptorSize);
	srvDesc.Format = waterTex->GetDesc().Format;
	md3dDevice->CreateShaderResourceView(waterTex.Get(), &srvDesc, hDescriptor);

	// next descriptor
	hDescriptor.Offset(1, mCbvSrvDescriptorSize);
	srvDesc.Format = fenceTex->GetDesc().Format;
	md3dDevice->CreateShaderResourceView(fenceTex.Get(), &srvDesc, hDescriptor);

	// next descriptor
	hDescriptor.Offset(1, mCbvSrvDescriptorSize);
	srvDesc.Format = brickTex->GetDesc().Format;
	md3dDevice->CreateShaderResourceView(brickTex.Get(), &srvDesc, hDescriptor);

	// next descriptor
	hDescriptor.Offset(1, mCbvSrvDescriptorSize);

	srvDesc.Format = brickDarkTex->GetDesc().Format;
	md3dDevice->CreateShaderResourceView(brickDarkTex.Get(), &srvDesc, hDescriptor);

	// next descriptor
	hDescriptor.Offset(1, mCbvSrvDescriptorSize);

	srvDesc.Format = bricks3->GetDesc().Format;
	md3dDevice->CreateShaderResourceView(bricks3.Get(), &srvDesc, hDescriptor);

	// next descriptor
	hDescriptor.Offset(1, mCbvSrvDescriptorSize);

	srvDesc.Format = stone->GetDesc().Format;
	md3dDevice->CreateShaderResourceView(stone.Get(), &srvDesc, hDescriptor);

	// next descriptor
	hDescriptor.Offset(1, mCbvSrvDescriptorSize);

	srvDesc.Format = black->GetDesc().Format;
	md3dDevice->CreateShaderResourceView(black.Get(), &srvDesc, hDescriptor);

	// next descriptor
	hDescriptor.Offset(1, mCbvSrvDescriptorSize);

	srvDesc.Format = sims->GetDesc().Format;
	md3dDevice->CreateShaderResourceView(sims.Get(), &srvDesc, hDescriptor);

	// next descriptor
	hDescriptor.Offset(1, mCbvSrvDescriptorSize);

	srvDesc.Format = bricksOnAcid->GetDesc().Format;
	md3dDevice->CreateShaderResourceView(bricksOnAcid.Get(), &srvDesc, hDescriptor);

	// next descriptor
	hDescriptor.Offset(1, mCbvSrvDescriptorSize);

	auto desc = treeArrayTex->GetDesc();
	srvDesc.ViewDimension = D3D12_SRV_DIMENSION_TEXTURE2DARRAY;
	srvDesc.Format = treeArrayTex->GetDesc().Format;
	srvDesc.Texture2DArray.MostDetailedMip = 0;
	srvDesc.Texture2DArray.MipLevels = -1;
	srvDesc.Texture2DArray.FirstArraySlice = 0;
	srvDesc.Texture2DArray.ArraySize = treeArrayTex->GetDesc().DepthOrArraySize;
	md3dDevice->CreateShaderResourceView(treeArrayTex.Get(), &srvDesc, hDescriptor);
}

void Castle::BuildShadersAndInputLayouts()
{
	const D3D_SHADER_MACRO defines[] =
	{
		NULL, NULL
	};

	const D3D_SHADER_MACRO alphaTestDefines[] =
	{
		"ALPHA_TEST", "1",
		NULL, NULL
	};

	mShaders["standardVS"] = d3dUtil::CompileShader(L"Shaders\\Default.hlsl", nullptr, "VS", "vs_5_0");
	mShaders["opaquePS"] = d3dUtil::CompileShader(L"Shaders\\Default.hlsl", defines, "PS", "ps_5_0");
	mShaders["alphaTestedPS"] = d3dUtil::CompileShader(L"Shaders\\Default.hlsl", alphaTestDefines, "PS", "ps_5_0");

	mShaders["treeSpriteVS"] = d3dUtil::CompileShader(L"Shaders\\TreeSprite.hlsl", nullptr, "VS", "vs_5_0");
	mShaders["treeSpriteGS"] = d3dUtil::CompileShader(L"Shaders\\TreeSprite.hlsl", nullptr, "GS", "gs_5_0");
	mShaders["treeSpritePS"] = d3dUtil::CompileShader(L"Shaders\\TreeSprite.hlsl", alphaTestDefines, "PS", "ps_5_0");

	mStdInputLayout =
	{
		{ "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
		{ "NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
		{ "TEXCOORD", 0, DXGI_FORMAT_R32G32_FLOAT, 0, 24, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
	};

	mTreeSpriteInputLayout =
	{
		{ "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
		{ "SIZE", 0, DXGI_FORMAT_R32G32_FLOAT, 0, 12, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
	};
}

void Castle::BuildLandGeometry()
{
	GeometryGenerator geoGen;
	GeometryGenerator::MeshData grid = geoGen.CreateGrid(320.0f, 160.0f, 50, 50, 80.0f, 0.0f);

	//
	// Extract the vertex elements we are interested and apply the height function to
	// each vertex.  In addition, color the vertices based on their height so we have
	// sandy looking beaches, grassy low hills, and snow mountain peaks.
	//

	std::vector<Vertex> vertices(grid.Vertices.size());
	for (size_t i = 0; i < grid.Vertices.size(); ++i)
	{
		auto& p = grid.Vertices[i].Position;
		vertices[i].Pos = p;
		//vertices[i].Pos.y = GetHillsHeight(p.x, p.z);
		float xPos = vertices[i].Pos.x;
		float zPos = vertices[i].Pos.z;
		if (xPos < -40 || xPos > 40 || zPos < -40 || zPos > 40)
		{
			vertices[i].Pos.y = 1.5f;
		}
		else if ((xPos > -30 && xPos < 30) && (zPos > -30 && zPos < 30))
		{
			vertices[i].Pos.y = 1.5f;
		}
		else
		{
			// ditch
			vertices[i].Pos.y = -20;
		}
		vertices[i].Normal = GetHillsNormal(p.x, p.z);
		vertices[i].TexC = grid.Vertices[i].TexC;
	}

	const UINT vbByteSize = (UINT)vertices.size() * sizeof(Vertex);

	std::vector<std::uint16_t> indices = grid.GetIndices16();
	const UINT ibByteSize = (UINT)indices.size() * sizeof(std::uint16_t);

	auto geo = std::make_unique<MeshGeometry>();
	geo->Name = "landGeo";

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

	SubmeshGeometry submesh;
	submesh.IndexCount = (UINT)indices.size();
	submesh.StartIndexLocation = 0;
	submesh.BaseVertexLocation = 0;

	geo->DrawArgs["grid"] = submesh;

	mGeometries["landGeo"] = std::move(geo);
}

void Castle::BuildWavesGeometry()
{
	std::vector<std::uint16_t> indices(3 * mWaves->TriangleCount()); // 3 indices per face
	assert(mWaves->VertexCount() < 0x0000ffff);

	// Iterate over each quad.
	int m = mWaves->RowCount();
	int n = mWaves->ColumnCount();
	int k = 0;
	for (int i = 0; i < m - 1; ++i)
	{
		for (int j = 0; j < n - 1; ++j)
		{
			indices[k] = i * n + j;
			indices[k + 1] = i * n + j + 1;
			indices[k + 2] = (i + 1)*n + j;

			indices[k + 3] = (i + 1)*n + j;
			indices[k + 4] = i * n + j + 1;
			indices[k + 5] = (i + 1)*n + j + 1;

			k += 6; // next quad
		}
	}

	UINT vbByteSize = mWaves->VertexCount() * sizeof(Vertex);
	UINT ibByteSize = (UINT)indices.size() * sizeof(std::uint16_t);

	auto geo = std::make_unique<MeshGeometry>();
	geo->Name = "waterGeo";

	// Set dynamically.
	geo->VertexBufferCPU = nullptr;
	geo->VertexBufferGPU = nullptr;

	ThrowIfFailed(D3DCreateBlob(ibByteSize, &geo->IndexBufferCPU));
	CopyMemory(geo->IndexBufferCPU->GetBufferPointer(), indices.data(), ibByteSize);

	geo->IndexBufferGPU = d3dUtil::CreateDefaultBuffer(md3dDevice.Get(),
		mCommandList.Get(), indices.data(), ibByteSize, geo->IndexBufferUploader);

	geo->VertexByteStride = sizeof(Vertex);
	geo->VertexBufferByteSize = vbByteSize;
	geo->IndexFormat = DXGI_FORMAT_R16_UINT;
	geo->IndexBufferByteSize = ibByteSize;

	SubmeshGeometry submesh;
	submesh.IndexCount = (UINT)indices.size();
	submesh.StartIndexLocation = 0;
	submesh.BaseVertexLocation = 0;

	geo->DrawArgs["grid"] = submesh;

	mGeometries["waterGeo"] = std::move(geo);
}

void Castle::BuildBoxGeometry()
{
	GeometryGenerator geoGen;
	GeometryGenerator::MeshData box = geoGen.CreateBox(8.0f, 8.0f, 8.0f, 3);

	std::vector<Vertex> vertices(box.Vertices.size());
	for (size_t i = 0; i < box.Vertices.size(); ++i)
	{
		auto& p = box.Vertices[i].Position;
		vertices[i].Pos = p;
		vertices[i].Normal = box.Vertices[i].Normal;
		vertices[i].TexC = box.Vertices[i].TexC;
	}

	const UINT vbByteSize = (UINT)vertices.size() * sizeof(Vertex);

	std::vector<std::uint16_t> indices = box.GetIndices16();
	const UINT ibByteSize = (UINT)indices.size() * sizeof(std::uint16_t);

	auto geo = std::make_unique<MeshGeometry>();
	geo->Name = "boxGeo";

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

	SubmeshGeometry submesh;
	submesh.IndexCount = (UINT)indices.size();
	submesh.StartIndexLocation = 0;
	submesh.BaseVertexLocation = 0;

	geo->DrawArgs["box"] = submesh;

	mGeometries["boxGeo"] = std::move(geo);
}

void Castle::BuildShapeGeometry()
{
	GeometryGenerator geoGen;
	GeometryGenerator::MeshData box = geoGen.CreateBox(1.5f, 1.5f, 1.5f, 3);
	GeometryGenerator::MeshData grid = geoGen.CreateGrid(20.0f, 30.0f, 60, 40, 0, 0);
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
		vertices[k].TexC = box.Vertices[i].TexC;
	}

	for (size_t i = 0; i < grid.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = grid.Vertices[i].Position;
		vertices[k].Normal = grid.Vertices[i].Normal;
		vertices[k].TexC = grid.Vertices[i].TexC;
	}

	for (size_t i = 0; i < sphere.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = sphere.Vertices[i].Position;
		vertices[k].Normal = sphere.Vertices[i].Normal;
		vertices[k].TexC = sphere.Vertices[i].TexC;
	}

	for (size_t i = 0; i < cylinder.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = cylinder.Vertices[i].Position;
		vertices[k].Normal = cylinder.Vertices[i].Normal;
		vertices[k].TexC = cylinder.Vertices[i].TexC;
	}

	for (size_t i = 0; i < star.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = star.Vertices[i].Position;
		vertices[k].TexC = star.Vertices[i].TexC;
		vertices[k].Normal = star.Vertices[i].Normal;
	}

	for (size_t i = 0; i < wedge.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = wedge.Vertices[i].Position;
		vertices[k].Normal = wedge.Vertices[i].Normal;
		vertices[k].TexC = wedge.Vertices[i].TexC;
	}

	for (size_t i = 0; i < triPrism.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = triPrism.Vertices[i].Position;
		vertices[k].Normal = triPrism.Vertices[i].Normal;
		vertices[k].TexC = triPrism.Vertices[i].TexC;
	}

	for (size_t i = 0; i < cube.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = cube.Vertices[i].Position;
		vertices[k].Normal = cube.Vertices[i].Normal;
		vertices[k].TexC = cube.Vertices[i].TexC;
	}

	for (size_t i = 0; i < pyramid.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = pyramid.Vertices[i].Position;
		vertices[k].Normal = pyramid.Vertices[i].Normal;
		vertices[k].TexC = pyramid.Vertices[i].TexC;
	}

	for (size_t i = 0; i < hexagon.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = hexagon.Vertices[i].Position;
		vertices[k].Normal = hexagon.Vertices[i].Normal;
		vertices[k].TexC = hexagon.Vertices[i].TexC;
	}

	for (size_t i = 0; i < simsInd.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = simsInd.Vertices[i].Position;
		vertices[k].Normal = simsInd.Vertices[i].Normal;
		vertices[k].TexC = simsInd.Vertices[i].TexC;
	}

	for (size_t i = 0; i < myCylinder.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = myCylinder.Vertices[i].Position;
		vertices[k].Normal = myCylinder.Vertices[i].Normal;
		vertices[k].TexC = myCylinder.Vertices[i].TexC;
	}

	for (size_t i = 0; i < cone.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = cone.Vertices[i].Position;
		vertices[k].Normal = cone.Vertices[i].Normal;
		vertices[k].TexC = cone.Vertices[i].TexC;
	}

	for (size_t i = 0; i < mySphere.Vertices.size(); ++i, ++k)
	{
		vertices[k].Pos = mySphere.Vertices[i].Position;
		vertices[k].Normal = mySphere.Vertices[i].Normal;
		vertices[k].TexC = mySphere.Vertices[i].TexC;
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

	/* Our shapes */
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

void Castle::BuildTreeSpritesGeometry()
{
	//step5
	struct TreeSpriteVertex
	{
		XMFLOAT3 Pos;
		XMFLOAT2 Size;
	};

	static const int treeCount = 20;
	std::array<TreeSpriteVertex, treeCount> vertices;

	vertices[0].Pos  = XMFLOAT3(-75, 10,  75);
	vertices[1].Pos  = XMFLOAT3(-45, 10,  75);
	vertices[2].Pos  = XMFLOAT3(-15, 10,  75);
	vertices[3].Pos  = XMFLOAT3( 15, 10,  75);
	vertices[4].Pos  = XMFLOAT3( 45, 10,  75);
	vertices[5].Pos  = XMFLOAT3( 75, 10,  75);
	vertices[6].Pos  = XMFLOAT3( 105, 10,  75);
	vertices[7].Pos  = XMFLOAT3( 145, 10,  75);

	vertices[8].Pos  = XMFLOAT3( 145, 10, -75);
	vertices[9].Pos  = XMFLOAT3( 105, 10, -75);
	vertices[10].Pos = XMFLOAT3( 75, 10, -75);
	vertices[11].Pos = XMFLOAT3( 45, 10, -75);
	vertices[12].Pos = XMFLOAT3( 15, 10, -75);
	vertices[13].Pos = XMFLOAT3(-15, 10, -75);
	vertices[14].Pos = XMFLOAT3(-45, 10, -75);
	vertices[15].Pos = XMFLOAT3(-75, 10, -75);
	vertices[16].Pos = XMFLOAT3(-75, 10, -45);
	vertices[17].Pos = XMFLOAT3(-75, 10, -15);
	vertices[18].Pos = XMFLOAT3(-75, 10,  15);
	vertices[19].Pos = XMFLOAT3(-75, 10,  45);

	for (UINT i = 0; i < treeCount; ++i)
	{
		vertices[i].Size = XMFLOAT2(20.0f, 20.0f);
	}

	std::array<std::uint16_t, treeCount> indices = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};

	const UINT vbByteSize = (UINT)vertices.size() * sizeof(TreeSpriteVertex);
	const UINT ibByteSize = (UINT)indices.size() * sizeof(std::uint16_t);

	auto geo = std::make_unique<MeshGeometry>();
	geo->Name = "treeSpritesGeo";

	ThrowIfFailed(D3DCreateBlob(vbByteSize, &geo->VertexBufferCPU));
	CopyMemory(geo->VertexBufferCPU->GetBufferPointer(), vertices.data(), vbByteSize);

	ThrowIfFailed(D3DCreateBlob(ibByteSize, &geo->IndexBufferCPU));
	CopyMemory(geo->IndexBufferCPU->GetBufferPointer(), indices.data(), ibByteSize);

	geo->VertexBufferGPU = d3dUtil::CreateDefaultBuffer(md3dDevice.Get(),
		mCommandList.Get(), vertices.data(), vbByteSize, geo->VertexBufferUploader);

	geo->IndexBufferGPU = d3dUtil::CreateDefaultBuffer(md3dDevice.Get(),
		mCommandList.Get(), indices.data(), ibByteSize, geo->IndexBufferUploader);

	geo->VertexByteStride = sizeof(TreeSpriteVertex);
	geo->VertexBufferByteSize = vbByteSize;
	geo->IndexFormat = DXGI_FORMAT_R16_UINT;
	geo->IndexBufferByteSize = ibByteSize;

	SubmeshGeometry submesh;
	submesh.IndexCount = (UINT)indices.size();
	submesh.StartIndexLocation = 0;
	submesh.BaseVertexLocation = 0;

	geo->DrawArgs["points"] = submesh;

	mGeometries["treeSpritesGeo"] = std::move(geo);
}

void Castle::BuildPSOs()
{
	D3D12_GRAPHICS_PIPELINE_STATE_DESC opaquePsoDesc;

	//
	// PSO for opaque objects.
	//
	ZeroMemory(&opaquePsoDesc, sizeof(D3D12_GRAPHICS_PIPELINE_STATE_DESC));
	opaquePsoDesc.InputLayout = { mStdInputLayout.data(), (UINT)mStdInputLayout.size() };
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

	//there is abug with F2 key that is supposed to turn on the multisampling!
//Set4xMsaaState(true);
	//m4xMsaaState = true;

	opaquePsoDesc.SampleDesc.Count = m4xMsaaState ? 4 : 1;
	opaquePsoDesc.SampleDesc.Quality = m4xMsaaState ? (m4xMsaaQuality - 1) : 0;
	opaquePsoDesc.DSVFormat = mDepthStencilFormat;
	ThrowIfFailed(md3dDevice->CreateGraphicsPipelineState(&opaquePsoDesc, IID_PPV_ARGS(&mPSOs["opaque"])));

	//
	// PSO for transparent objects
	//

	D3D12_GRAPHICS_PIPELINE_STATE_DESC transparentPsoDesc = opaquePsoDesc;

	D3D12_RENDER_TARGET_BLEND_DESC transparencyBlendDesc;
	transparencyBlendDesc.BlendEnable = true;
	transparencyBlendDesc.LogicOpEnable = false;
	transparencyBlendDesc.SrcBlend = D3D12_BLEND_SRC_ALPHA;
	transparencyBlendDesc.DestBlend = D3D12_BLEND_INV_SRC_ALPHA;
	transparencyBlendDesc.BlendOp = D3D12_BLEND_OP_ADD;
	transparencyBlendDesc.SrcBlendAlpha = D3D12_BLEND_ONE;
	transparencyBlendDesc.DestBlendAlpha = D3D12_BLEND_ZERO;
	transparencyBlendDesc.BlendOpAlpha = D3D12_BLEND_OP_ADD;
	transparencyBlendDesc.LogicOp = D3D12_LOGIC_OP_NOOP;
	transparencyBlendDesc.RenderTargetWriteMask = D3D12_COLOR_WRITE_ENABLE_ALL;

	//transparentPsoDesc.BlendState.AlphaToCoverageEnable = true;

	transparentPsoDesc.BlendState.RenderTarget[0] = transparencyBlendDesc;
	ThrowIfFailed(md3dDevice->CreateGraphicsPipelineState(&transparentPsoDesc, IID_PPV_ARGS(&mPSOs["transparent"])));

	//
	// PSO for alpha tested objects
	//

	D3D12_GRAPHICS_PIPELINE_STATE_DESC alphaTestedPsoDesc = opaquePsoDesc;
	alphaTestedPsoDesc.PS =
	{
		reinterpret_cast<BYTE*>(mShaders["alphaTestedPS"]->GetBufferPointer()),
		mShaders["alphaTestedPS"]->GetBufferSize()
	};
	alphaTestedPsoDesc.RasterizerState.CullMode = D3D12_CULL_MODE_NONE;
	ThrowIfFailed(md3dDevice->CreateGraphicsPipelineState(&alphaTestedPsoDesc, IID_PPV_ARGS(&mPSOs["alphaTested"])));

	//
	// PSO for tree sprites
	//
	D3D12_GRAPHICS_PIPELINE_STATE_DESC treeSpritePsoDesc = opaquePsoDesc;
	treeSpritePsoDesc.VS =
	{
		reinterpret_cast<BYTE*>(mShaders["treeSpriteVS"]->GetBufferPointer()),
		mShaders["treeSpriteVS"]->GetBufferSize()
	};
	treeSpritePsoDesc.GS =
	{
		reinterpret_cast<BYTE*>(mShaders["treeSpriteGS"]->GetBufferPointer()),
		mShaders["treeSpriteGS"]->GetBufferSize()
	};
	treeSpritePsoDesc.PS =
	{
		reinterpret_cast<BYTE*>(mShaders["treeSpritePS"]->GetBufferPointer()),
		mShaders["treeSpritePS"]->GetBufferSize()
	};
	//step1
	treeSpritePsoDesc.PrimitiveTopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_POINT;
	treeSpritePsoDesc.InputLayout = { mTreeSpriteInputLayout.data(), (UINT)mTreeSpriteInputLayout.size() };
	treeSpritePsoDesc.RasterizerState.CullMode = D3D12_CULL_MODE_NONE;

	ThrowIfFailed(md3dDevice->CreateGraphicsPipelineState(&treeSpritePsoDesc, IID_PPV_ARGS(&mPSOs["treeSprites"])));
}

void Castle::BuildFrameResources()
{
	for (int i = 0; i < gNumFrameResources; ++i)
	{
		mFrameResources.push_back(std::make_unique<FrameResource>(md3dDevice.Get(),
			1, (UINT)mAllRitems.size(), (UINT)mMaterials.size(), mWaves->VertexCount()));
	}
}

void Castle::BuildRenderGeoItem(int index, std::string itemName, std::string material, XMFLOAT3 scaling, XMFLOAT3 rotation, XMFLOAT3 translation) {

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
	mRitemLayer[(int)RenderLayer::Opaque].push_back(genItem.get());
	mAllRitems.push_back(std::move(genItem));

}

void Castle::BuildMaterials()
{
	auto grass = std::make_unique<Material>();
	grass->Name = "grass";
	grass->MatCBIndex = 0;
	grass->DiffuseSrvHeapIndex = 0;
	grass->DiffuseAlbedo = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);
	grass->FresnelR0 = XMFLOAT3(0.01f, 0.01f, 0.01f);
	grass->Roughness = 0.125f;

	// This is not a good water material definition, but we do not have all the rendering
	// tools we need (transparency, environment reflection), so we fake it for now.
	auto water = std::make_unique<Material>();
	water->Name = "water";
	water->MatCBIndex = 1;
	water->DiffuseSrvHeapIndex = 1;
	water->DiffuseAlbedo = XMFLOAT4(1.0f, 1.0f, 1.0f, 0.5f);
	water->FresnelR0 = XMFLOAT3(0.1f, 0.1f, 0.1f);
	water->Roughness = 0.0f;

	auto wirefence = std::make_unique<Material>();
	wirefence->Name = "wirefence";
	wirefence->MatCBIndex = 2;
	wirefence->DiffuseSrvHeapIndex = 2;
	wirefence->DiffuseAlbedo = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);
	wirefence->FresnelR0 = XMFLOAT3(0.02f, 0.02f, 0.02f);
	wirefence->Roughness = 0.25f;

	auto treeSprites = std::make_unique<Material>();
	treeSprites->Name = "treeSprites";
	treeSprites->MatCBIndex = 3;
	treeSprites->DiffuseSrvHeapIndex = 10;
	treeSprites->DiffuseAlbedo = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);
	treeSprites->FresnelR0 = XMFLOAT3(0.01f, 0.01f, 0.01f);
	treeSprites->Roughness = 0.125f;

	auto brick = std::make_unique<Material>();
	brick->Name = "brick";
	brick->MatCBIndex = 4;
	brick->DiffuseSrvHeapIndex = 3;
	brick->DiffuseAlbedo = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);
	brick->FresnelR0 = XMFLOAT3(0.01f, 0.01f, 0.01f);
	brick->Roughness = 0.125f;

	auto brickDark = std::make_unique<Material>();
	brickDark->Name = "brickDark";
	brickDark->MatCBIndex = 5;
	brickDark->DiffuseSrvHeapIndex = 4;
	brickDark->DiffuseAlbedo = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);
	brickDark->FresnelR0 = XMFLOAT3(0.01f, 0.01f, 0.01f);
	brickDark->Roughness = 0.125f;

	auto bricks3 = std::make_unique<Material>();
	bricks3->Name = "bricks3";
	bricks3->MatCBIndex = 6;
	bricks3->DiffuseSrvHeapIndex = 5;
	bricks3->DiffuseAlbedo = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);
	bricks3->FresnelR0 = XMFLOAT3(0.01f, 0.01f, 0.01f);
	bricks3->Roughness = 0.125f;

	auto stone = std::make_unique<Material>();
	stone->Name = "stone";
	stone->MatCBIndex = 7;
	stone->DiffuseSrvHeapIndex = 6;
	stone->DiffuseAlbedo = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);
	stone->FresnelR0 = XMFLOAT3(0.01f, 0.01f, 0.01f);
	stone->Roughness = 0.125f;

	auto black = std::make_unique<Material>();
	black->Name = "black";
	black->MatCBIndex = 8;
	black->DiffuseSrvHeapIndex = 7;
	black->DiffuseAlbedo = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);
	black->FresnelR0 = XMFLOAT3(0.01f, 0.01f, 0.01f);
	black->Roughness = 0.125f;

	auto sims = std::make_unique<Material>();
	sims->Name = "sims";
	sims->MatCBIndex = 9;
	sims->DiffuseSrvHeapIndex = 8;
	sims->DiffuseAlbedo = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);
	sims->FresnelR0 = XMFLOAT3(0.01f, 0.01f, 0.01f);
	sims->Roughness = 0.125f;

	auto bricksOnAcid = std::make_unique<Material>();
	bricksOnAcid->Name = "bricksOnAcid";
	bricksOnAcid->MatCBIndex = 10;
	bricksOnAcid->DiffuseSrvHeapIndex = 9;
	bricksOnAcid->DiffuseAlbedo = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);
	bricksOnAcid->FresnelR0 = XMFLOAT3(0.01f, 0.01f, 0.01f);
	bricksOnAcid->Roughness = 0.125f;

	mMaterials["grass"] = std::move(grass);
	mMaterials["water"] = std::move(water);
	mMaterials["wirefence"] = std::move(wirefence);
	mMaterials["treeSprites"] = std::move(treeSprites);
	mMaterials["brick"] = std::move(brick);
	mMaterials["brickDark"] = std::move(brickDark);
	mMaterials["bricks3"] = std::move(bricks3);
	mMaterials["stone"] = std::move(stone);
	mMaterials["black"] = std::move(black);
	mMaterials["sims"] = std::move(sims);
	mMaterials["bricksOnAcid"] = std::move(bricksOnAcid);
}

void Castle::BuildRenderGeoItems()
{
	int i = 0;
	auto wavesRitem = std::make_unique<RenderItem>();
	XMMATRIX transResult = XMMatrixScaling(5.0f, 5.0f, 1.0f);
	transResult *= XMMatrixTranslation(0, -1.5f, 0);
	//XMStoreFloat4x4(&wavesRitem->World, transResult);

	XMStoreFloat4x4(&wavesRitem->TexTransform, transResult);
	wavesRitem->ObjCBIndex = i++;
	wavesRitem->Mat = mMaterials["water"].get();
	wavesRitem->Geo = mGeometries["waterGeo"].get();
	wavesRitem->PrimitiveType = D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST;
	wavesRitem->IndexCount = wavesRitem->Geo->DrawArgs["grid"].IndexCount;
	wavesRitem->StartIndexLocation = wavesRitem->Geo->DrawArgs["grid"].StartIndexLocation;
	wavesRitem->BaseVertexLocation = wavesRitem->Geo->DrawArgs["grid"].BaseVertexLocation;
	mWavesRitem = wavesRitem.get();
	mRitemLayer[(int)RenderLayer::Transparent].push_back(wavesRitem.get());
	mAllRitems.push_back(std::move(wavesRitem));

	auto gridRitem = std::make_unique<RenderItem>();
	gridRitem->World = MathHelper::Identity4x4();
	XMStoreFloat4x4(&gridRitem->TexTransform, XMMatrixScaling(5.0f, 5.0f, 1.0f));
	gridRitem->ObjCBIndex = i++;
	gridRitem->Mat = mMaterials["grass"].get();
	gridRitem->Geo = mGeometries["landGeo"].get();
	gridRitem->PrimitiveType = D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST;
	gridRitem->IndexCount = gridRitem->Geo->DrawArgs["grid"].IndexCount;
	gridRitem->StartIndexLocation = gridRitem->Geo->DrawArgs["grid"].StartIndexLocation;
	gridRitem->BaseVertexLocation = gridRitem->Geo->DrawArgs["grid"].BaseVertexLocation;
	mRitemLayer[(int)RenderLayer::Opaque].push_back(gridRitem.get());
	mAllRitems.push_back(std::move(gridRitem));

	auto treeSpritesRitem = std::make_unique<RenderItem>();
	treeSpritesRitem->World = MathHelper::Identity4x4();
	treeSpritesRitem->ObjCBIndex = i++;
	treeSpritesRitem->Mat = mMaterials["treeSprites"].get();
	treeSpritesRitem->Geo = mGeometries["treeSpritesGeo"].get();
	//step2
	treeSpritesRitem->PrimitiveType = D3D_PRIMITIVE_TOPOLOGY_POINTLIST;
	treeSpritesRitem->IndexCount = treeSpritesRitem->Geo->DrawArgs["points"].IndexCount;
	treeSpritesRitem->StartIndexLocation = treeSpritesRitem->Geo->DrawArgs["points"].StartIndexLocation;
	treeSpritesRitem->BaseVertexLocation = treeSpritesRitem->Geo->DrawArgs["points"].BaseVertexLocation;

	mRitemLayer[(int)RenderLayer::AlphaTestedTreeSprites].push_back(treeSpritesRitem.get());
	mAllRitems.push_back(std::move(treeSpritesRitem));

	// The castle
	//// Bridge
	BuildRenderGeoItem(i++, "cube", "stone", XMFLOAT3(20.05f, 1.5, 8.33f), XMFLOAT3(0, 0, 0), XMFLOAT3(30.74f, 1, 0));
	//// Walls
	BuildRenderGeoItem(i++, "cube", "brick", XMFLOAT3(27.86f, 1.76f, 44.744f), XMFLOAT3(0, 180, 90), XMFLOAT3(-22.27f, 2, 0.144f));
	BuildRenderGeoItem(i++, "cube", "brick", XMFLOAT3(27.86f, 1.76f, 44.744f), XMFLOAT3(0, 0, 90), XMFLOAT3(20.74f, 2, 0.144f));
	BuildRenderGeoItem(i++, "cube", "brick", XMFLOAT3(27.86f, 1.76f, 44.744f), XMFLOAT3(0, 90, 90), XMFLOAT3(-0.72f, 2, -21.68f));
	BuildRenderGeoItem(i++, "cube", "brick", XMFLOAT3(27.86f, 1.76f, 44.744f), XMFLOAT3(0, 270, 90), XMFLOAT3(-0.72f, 2, 21.49f));
	//// Central Cube
	BuildRenderGeoItem(i++, "cube", "brick", XMFLOAT3(23, 15, 23), XMFLOAT3(0, 180, 90), XMFLOAT3(0, 4.4f, 0));
	//// Portal
	BuildRenderGeoItem(i++, "cube", "black", XMFLOAT3(19.05f, 1.95f, 8.32f), XMFLOAT3(0, 0, 0), XMFLOAT3(20.76f, -2.4f, 0));
	BuildRenderGeoItem(i++, "myCylinder", "black", XMFLOAT3(8.38f, 0.96f, 8.38f), XMFLOAT3(0, 0, 0), XMFLOAT3(20.76f, 7.02f, 0.011f));
	//// Towers
	BuildRenderGeoItem(i++, "myCylinder", "brickDark", XMFLOAT3(4, 40, 4), XMFLOAT3(0, 0, 0), XMFLOAT3(20.8f, 2.73f, -21.5f));
	BuildRenderGeoItem(i++, "myCylinder", "brickDark", XMFLOAT3(4, 40, 4), XMFLOAT3(0, 0, 0), XMFLOAT3(-22.21f, 2.73f, -21.5f));
	BuildRenderGeoItem(i++, "myCylinder", "brickDark", XMFLOAT3(4, 40, 4), XMFLOAT3(0, 0, 0), XMFLOAT3(20.8f, 2.73f, 21.67f));
	BuildRenderGeoItem(i++, "myCylinder", "brickDark", XMFLOAT3(4, 40, 4), XMFLOAT3(0, 0, 0), XMFLOAT3(-20.19f, 2.73f, 21.67f));
	BuildRenderGeoItem(i++, "myCylinder", "brickDark", XMFLOAT3(5, 0.3f, 5), XMFLOAT3(0, 0, 0), XMFLOAT3(20.8f, 22.83f, -21.5f));
	BuildRenderGeoItem(i++, "myCylinder", "brickDark", XMFLOAT3(5, 0.3f, 5), XMFLOAT3(0, 0, 0), XMFLOAT3(-22.21f, 22.83f, -21.5f));
	BuildRenderGeoItem(i++, "myCylinder", "brickDark", XMFLOAT3(5, 0.3f, 5), XMFLOAT3(0, 0, 0), XMFLOAT3(20.8f, 22.83f, 21.67f));
	BuildRenderGeoItem(i++, "myCylinder", "brickDark", XMFLOAT3(5, 0.3f, 5), XMFLOAT3(0, 0, 0), XMFLOAT3(-20.19f, 22.83f, 21.67f));
	//// Strange sphere on top of central building
	BuildRenderGeoItem(i++, "mySphere", "brickDark", XMFLOAT3(15, 15, 15), XMFLOAT3(0, 0, 0), XMFLOAT3(0, 11.91f, 0));
	//// Strange indicator
	BuildRenderGeoItem(i++, "simsInd", "sims", XMFLOAT3(5, 1, 1), XMFLOAT3(90, 0, 0), XMFLOAT3(21.55f, 5.9f, -11));
	BuildRenderGeoItem(i++, "simsInd", "sims", XMFLOAT3(5, 1, 1), XMFLOAT3(90, 0, 0), XMFLOAT3(21.55f, 5.9f, 11));
	//// StarPortal
	BuildRenderGeoItem(i++, "star", "black", XMFLOAT3(15, 15, 1), XMFLOAT3(0, 90, 0), XMFLOAT3(21.55f, 5.9f, 0));
	//// Pyramids for good luck
	BuildRenderGeoItem(i++, "pyramid", "bricks3", XMFLOAT3(3, 5, 3), XMFLOAT3(0, 0, 0), XMFLOAT3(20.8f, 24.83f, -21.5f));
	BuildRenderGeoItem(i++, "pyramid", "bricks3", XMFLOAT3(3, 5, 3), XMFLOAT3(0, 0, 0), XMFLOAT3(-22.21f, 24.83f, -21.5f));
	BuildRenderGeoItem(i++, "pyramid", "bricks3", XMFLOAT3(3, 5, 3), XMFLOAT3(0, 0, 0), XMFLOAT3(20.8f, 24.83f, 21.67f));
	BuildRenderGeoItem(i++, "pyramid", "bricks3", XMFLOAT3(3, 5, 3), XMFLOAT3(0, 0, 0), XMFLOAT3(-20.19f, 24.83f, 21.67f));
	//// The spaceship
	BuildRenderGeoItem(i++, "hexagon", "bricksOnAcid", XMFLOAT3(50, 5, 50), XMFLOAT3(0, 0, 0), XMFLOAT3(0, 50, 0));
	//// Pathways

	BuildRenderGeoItem(i++, "wedge", "bricks3", XMFLOAT3(2, 2, 2), XMFLOAT3(0,  90, 0), XMFLOAT3(mBasePathway    , 2.26f, mPathwayLeft));
	BuildRenderGeoItem(i++, "wedge", "bricks3", XMFLOAT3(2, 2, 2), XMFLOAT3(0,  90, 0), XMFLOAT3(mBasePathway + 2, 2.26f, mPathwayLeft));
	BuildRenderGeoItem(i++, "wedge", "bricks3", XMFLOAT3(2, 2, 2), XMFLOAT3(0,  90, 0), XMFLOAT3(mBasePathway + 4, 2.26f, mPathwayLeft));
	BuildRenderGeoItem(i++, "wedge", "bricks3", XMFLOAT3(2, 2, 2), XMFLOAT3(0,  90, 0), XMFLOAT3(mBasePathway + 6, 2.26f, mPathwayLeft));
	BuildRenderGeoItem(i++, "wedge", "bricks3", XMFLOAT3(2, 2, 2), XMFLOAT3(0,  90, 0), XMFLOAT3(mBasePathway + 8, 2.26f, mPathwayLeft));
	BuildRenderGeoItem(i++, "wedge", "bricks3", XMFLOAT3(2, 2, 2), XMFLOAT3(0, -90, 0), XMFLOAT3(mBasePathway    , 2.26f, mPathwayRight));
	BuildRenderGeoItem(i++, "wedge", "bricks3", XMFLOAT3(2, 2, 2), XMFLOAT3(0, -90, 0), XMFLOAT3(mBasePathway + 2, 2.26f, mPathwayRight));
	BuildRenderGeoItem(i++, "wedge", "bricks3", XMFLOAT3(2, 2, 2), XMFLOAT3(0, -90, 0), XMFLOAT3(mBasePathway + 4, 2.26f, mPathwayRight));
	BuildRenderGeoItem(i++, "wedge", "bricks3", XMFLOAT3(2, 2, 2), XMFLOAT3(0, -90, 0), XMFLOAT3(mBasePathway + 6, 2.26f, mPathwayRight));
	BuildRenderGeoItem(i++, "wedge", "bricks3", XMFLOAT3(2, 2, 2), XMFLOAT3(0, -90, 0), XMFLOAT3(mBasePathway + 8, 2.26f, mPathwayRight));
	// Beacons
	BuildRenderGeoItem(i++, "myCylinder", "stone", XMFLOAT3(2, 7.5f, 2), XMFLOAT3(0, 0, 0), XMFLOAT3(mBasePathway + 10, 2.26f, mPathwayLeft + 1));
	BuildRenderGeoItem(i++, "myCylinder", "stone", XMFLOAT3(2, 7.5f, 2), XMFLOAT3(0, 0, 0), XMFLOAT3(mBasePathway + 10, 2.26f, mPathwayRight - 1));

	//// Why a prism here?
	BuildRenderGeoItem(i++, "triPrism", "bricks3", XMFLOAT3(15, 15, 15), XMFLOAT3(0, 0, 0), XMFLOAT3(0, 20, -3.0f));

	// The Maze
	BuildRenderGeoItem(i++, "cube", "bricks3", XMFLOAT3(30, 30, 30), XMFLOAT3(0, 0, 0), XMFLOAT3(100, 15 + 1.5f, 0));
	blocks[0] = BlockInfo(XMFLOAT3(100, 15 + 1.5f, 0), XMFLOAT3(30, 30, 30));
}

void Castle::DrawRenderItems(ID3D12GraphicsCommandList* cmdList, const std::vector<RenderItem*>& ritems)
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
		//step3
		cmdList->IASetPrimitiveTopology(ri->PrimitiveType);

		CD3DX12_GPU_DESCRIPTOR_HANDLE tex(mSrvDescriptorHeap->GetGPUDescriptorHandleForHeapStart());
		tex.Offset(ri->Mat->DiffuseSrvHeapIndex, mCbvSrvDescriptorSize);

		D3D12_GPU_VIRTUAL_ADDRESS objCBAddress = objectCB->GetGPUVirtualAddress() + ri->ObjCBIndex*objCBByteSize;
		D3D12_GPU_VIRTUAL_ADDRESS matCBAddress = matCB->GetGPUVirtualAddress() + ri->Mat->MatCBIndex*matCBByteSize;

		cmdList->SetGraphicsRootDescriptorTable(0, tex);
		cmdList->SetGraphicsRootConstantBufferView(1, objCBAddress);
		cmdList->SetGraphicsRootConstantBufferView(3, matCBAddress);

		cmdList->DrawIndexedInstanced(ri->IndexCount, 1, ri->StartIndexLocation, ri->BaseVertexLocation, 0);
	}
}

std::array<const CD3DX12_STATIC_SAMPLER_DESC, 6> Castle::GetStaticSamplers()
{
	// Applications usually only need a handful of samplers.  So just define them all up front
	// and keep them available as part of the root signature.  

	const CD3DX12_STATIC_SAMPLER_DESC pointWrap(
		0, // shaderRegister
		D3D12_FILTER_MIN_MAG_MIP_POINT, // filter
		D3D12_TEXTURE_ADDRESS_MODE_WRAP,  // addressU
		D3D12_TEXTURE_ADDRESS_MODE_WRAP,  // addressV
		D3D12_TEXTURE_ADDRESS_MODE_WRAP); // addressW

	const CD3DX12_STATIC_SAMPLER_DESC pointClamp(
		1, // shaderRegister
		D3D12_FILTER_MIN_MAG_MIP_POINT, // filter
		D3D12_TEXTURE_ADDRESS_MODE_CLAMP,  // addressU
		D3D12_TEXTURE_ADDRESS_MODE_CLAMP,  // addressV
		D3D12_TEXTURE_ADDRESS_MODE_CLAMP); // addressW

	const CD3DX12_STATIC_SAMPLER_DESC linearWrap(
		2, // shaderRegister
		D3D12_FILTER_MIN_MAG_MIP_LINEAR, // filter
		D3D12_TEXTURE_ADDRESS_MODE_WRAP,  // addressU
		D3D12_TEXTURE_ADDRESS_MODE_WRAP,  // addressV
		D3D12_TEXTURE_ADDRESS_MODE_WRAP); // addressW

	const CD3DX12_STATIC_SAMPLER_DESC linearClamp(
		3, // shaderRegister
		D3D12_FILTER_MIN_MAG_MIP_LINEAR, // filter
		D3D12_TEXTURE_ADDRESS_MODE_CLAMP,  // addressU
		D3D12_TEXTURE_ADDRESS_MODE_CLAMP,  // addressV
		D3D12_TEXTURE_ADDRESS_MODE_CLAMP); // addressW

	const CD3DX12_STATIC_SAMPLER_DESC anisotropicWrap(
		4, // shaderRegister
		D3D12_FILTER_ANISOTROPIC, // filter
		D3D12_TEXTURE_ADDRESS_MODE_WRAP,  // addressU
		D3D12_TEXTURE_ADDRESS_MODE_WRAP,  // addressV
		D3D12_TEXTURE_ADDRESS_MODE_WRAP,  // addressW
		0.0f,                             // mipLODBias
		8);                               // maxAnisotropy

	const CD3DX12_STATIC_SAMPLER_DESC anisotropicClamp(
		5, // shaderRegister
		D3D12_FILTER_ANISOTROPIC, // filter
		D3D12_TEXTURE_ADDRESS_MODE_CLAMP,  // addressU
		D3D12_TEXTURE_ADDRESS_MODE_CLAMP,  // addressV
		D3D12_TEXTURE_ADDRESS_MODE_CLAMP,  // addressW
		0.0f,                              // mipLODBias
		8);                                // maxAnisotropy

	return {
		pointWrap, pointClamp,
		linearWrap, linearClamp,
		anisotropicWrap, anisotropicClamp };
}

float Castle::GetHillsHeight(float x, float z)const
{
	return 0.3f*(z*sinf(0.1f*x) + x * cosf(0.1f*z));
}

XMFLOAT3 Castle::GetHillsNormal(float x, float z)const
{
	// n = (-df/dx, 1, -df/dz)
	XMFLOAT3 n(
		-0.03f*z*cosf(0.1f*x) - 0.3f*cosf(0.1f*z),
		1.0f,
		-0.3f*sinf(0.1f*x) + 0.03f*x*sinf(0.1f*z));

	XMVECTOR unitNormal = XMVector3Normalize(XMLoadFloat3(&n));
	XMStoreFloat3(&n, unitNormal);

	return n;
}
