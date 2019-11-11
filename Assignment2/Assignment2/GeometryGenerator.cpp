//***************************************************************************************
// GeometryGenerator.cpp by Frank Luna (C) 2011 All Rights Reserved.
//***************************************************************************************

#include "GeometryGenerator.h"
#include <algorithm>
#include <DirectXMath.h>
#define PI 3.14159265

using namespace DirectX;

float GeometryGenerator::EpsilonCorrector(float value) {
	return (value < 3 * FLT_EPSILON && value > -3 * FLT_EPSILON) ? 0.0f : value;
}

DirectX::XMFLOAT3 GeometryGenerator::CrossProduct(DirectX::XMFLOAT3 vect_A, DirectX::XMFLOAT3 vect_B)
{
	XMFLOAT3 cross_P;
	cross_P.x = vect_A.y * vect_B.z - vect_A.z * vect_B.y;
	cross_P.y = vect_A.x * vect_B.z - vect_A.z * vect_B.x;
	cross_P.z = vect_A.x * vect_B.y - vect_A.y * vect_B.x;
	float length = sqrt(cross_P.x * cross_P.x + cross_P.y * cross_P.y + cross_P.z * cross_P.z);
	cross_P.x = epsilonCorrector(cross_P.x / length);
	cross_P.y = epsilonCorrector(cross_P.y / length);
	cross_P.z = epsilonCorrector(cross_P.z / length);
	return cross_P;
}
////
//XMFLOAT3 GeometryGenerator::vertexNormalGenerator(XMFLOAT3 left, XMFLOAT3 center, XMFLOAT3 right) {
//	XMFLOAT3 A, B, Normal;
//	A.x = left.x - center.x;
//	A.y = left.y - center.y;
//	A.z = left.z - center.z;
//	B.x = right.x - center.x;
//	B.y = right.y - center.y;
//	B.z = right.z - center.z;
//	Normal = crossProduct(A, B);
//	return Normal;
//}


GeometryGenerator::MeshData GeometryGenerator::CreateBox(float width, float height, float depth, uint32 numSubdivisions)
{
    MeshData meshData;

    //
	// Create the vertices.
	//

	Vertex v[24];

	float w2 = 0.5f*width;
	float h2 = 0.5f*height;
	float d2 = 0.5f*depth;
    
	// Fill in the front face vertex data.
	v[0] = Vertex(-w2, -h2, -d2);
	v[1] = Vertex(-w2, +h2, -d2);
	v[2] = Vertex(+w2, +h2, -d2);
	v[3] = Vertex(+w2, -h2, -d2);

	// Fill in the back face vertex data.
	v[4] = Vertex(-w2, -h2, +d2);
	v[5] = Vertex(+w2, -h2, +d2);
	v[6] = Vertex(+w2, +h2, +d2);
	v[7] = Vertex(-w2, +h2, +d2);

	// Fill in the top face vertex data.
	v[8]  = Vertex(-w2, +h2, -d2);
	v[9]  = Vertex(-w2, +h2, +d2);
	v[10] = Vertex(+w2, +h2, +d2);
	v[11] = Vertex(+w2, +h2, -d2);

	// Fill in the bottom face vertex data.
	v[12] = Vertex(-w2, -h2, -d2);
	v[13] = Vertex(+w2, -h2, -d2);
	v[14] = Vertex(+w2, -h2, +d2);
	v[15] = Vertex(-w2, -h2, +d2);

	// Fill in the left face vertex data.
	v[16] = Vertex(-w2, -h2, +d2);
	v[17] = Vertex(-w2, +h2, +d2);
	v[18] = Vertex(-w2, +h2, -d2);
	v[19] = Vertex(-w2, -h2, -d2);

	// Fill in the right face vertex data.
	v[20] = Vertex(+w2, -h2, -d2);
	v[21] = Vertex(+w2, +h2, -d2);
	v[22] = Vertex(+w2, +h2, +d2);
	v[23] = Vertex(+w2, -h2, +d2);

	meshData.Vertices.assign(&v[0], &v[24]);
 
	//
	// Create the indices.
	//

	uint32 i[36];

	// Fill in the front face index data
	i[0] = 0; i[1] = 1; i[2] = 2;
	i[3] = 0; i[4] = 2; i[5] = 3;

	// Fill in the back face index data
	i[6] = 4; i[7]  = 5; i[8]  = 6;
	i[9] = 4; i[10] = 6; i[11] = 7;

	// Fill in the top face index data
	i[12] = 8; i[13] =  9; i[14] = 10;
	i[15] = 8; i[16] = 10; i[17] = 11;

	// Fill in the bottom face index data
	i[18] = 12; i[19] = 13; i[20] = 14;
	i[21] = 12; i[22] = 14; i[23] = 15;

	// Fill in the left face index data
	i[24] = 16; i[25] = 17; i[26] = 18;
	i[27] = 16; i[28] = 18; i[29] = 19;

	// Fill in the right face index data
	i[30] = 20; i[31] = 21; i[32] = 22;
	i[33] = 20; i[34] = 22; i[35] = 23;

	meshData.Indices32.assign(&i[0], &i[36]);

    // Put a cap on the number of subdivisions.
    numSubdivisions = std::min<uint32>(numSubdivisions, 6u);

    for(uint32 i = 0; i < numSubdivisions; ++i)
        Subdivide(meshData);

    return meshData;
}

GeometryGenerator::MeshData GeometryGenerator::CreateSphere(float radius, uint32 sliceCount, uint32 stackCount)
{
    MeshData meshData;

	//
	// Compute the vertices stating at the top pole and moving down the stacks.
	//

	// Poles: note that there will be texture coordinate distortion as there is
	// not a unique point on the texture map to assign to the pole when mapping
	// a rectangular texture onto a sphere.
	Vertex topVertex(0.0f, +radius, 0.0f, 0.0f, +1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f);
	Vertex bottomVertex(0.0f, -radius, 0.0f, 0.0f, -1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f);

	meshData.Vertices.push_back( topVertex );

	float phiStep   = XM_PI/stackCount;
	float thetaStep = 2.0f*XM_PI/sliceCount;

	// Compute vertices for each stack ring (do not count the poles as rings).
	for(uint32 i = 1; i <= stackCount-1; ++i)
	{
		float phi = i*phiStep;

		// Vertices of ring.
        for(uint32 j = 0; j <= sliceCount; ++j)
		{
			float theta = j*thetaStep;

			Vertex v;

			// spherical to cartesian
			v.Position.x = radius*sinf(phi)*cosf(theta);
			v.Position.y = radius*cosf(phi);
			v.Position.z = radius*sinf(phi)*sinf(theta);

			// Partial derivative of P with respect to theta
			v.TangentU.x = -radius*sinf(phi)*sinf(theta);
			v.TangentU.y = 0.0f;
			v.TangentU.z = +radius*sinf(phi)*cosf(theta);

			XMVECTOR T = XMLoadFloat3(&v.TangentU);
			XMStoreFloat3(&v.TangentU, XMVector3Normalize(T));

			XMVECTOR p = XMLoadFloat3(&v.Position);
			XMStoreFloat3(&v.Normal, XMVector3Normalize(p));

			v.TexC.x = theta / XM_2PI;
			v.TexC.y = phi / XM_PI;

			meshData.Vertices.push_back( v );
		}
	}

	meshData.Vertices.push_back( bottomVertex );

	//
	// Compute indices for top stack.  The top stack was written first to the vertex buffer
	// and connects the top pole to the first ring.
	//

    for(uint32 i = 1; i <= sliceCount; ++i)
	{
		meshData.Indices32.push_back(0);
		meshData.Indices32.push_back(i+1);
		meshData.Indices32.push_back(i);
	}
	
	//
	// Compute indices for inner stacks (not connected to poles).
	//

	// Offset the indices to the index of the first vertex in the first ring.
	// This is just skipping the top pole vertex.
    uint32 baseIndex = 1;
    uint32 ringVertexCount = sliceCount + 1;
	for(uint32 i = 0; i < stackCount-2; ++i)
	{
		for(uint32 j = 0; j < sliceCount; ++j)
		{
			meshData.Indices32.push_back(baseIndex + i*ringVertexCount + j);
			meshData.Indices32.push_back(baseIndex + i*ringVertexCount + j+1);
			meshData.Indices32.push_back(baseIndex + (i+1)*ringVertexCount + j);

			meshData.Indices32.push_back(baseIndex + (i+1)*ringVertexCount + j);
			meshData.Indices32.push_back(baseIndex + i*ringVertexCount + j+1);
			meshData.Indices32.push_back(baseIndex + (i+1)*ringVertexCount + j+1);
		}
	}

	//
	// Compute indices for bottom stack.  The bottom stack was written last to the vertex buffer
	// and connects the bottom pole to the bottom ring.
	//

	// South pole vertex was added last.
	uint32 southPoleIndex = (uint32)meshData.Vertices.size()-1;

	// Offset the indices to the index of the first vertex in the last ring.
	baseIndex = southPoleIndex - ringVertexCount;
	
	for(uint32 i = 0; i < sliceCount; ++i)
	{
		meshData.Indices32.push_back(southPoleIndex);
		meshData.Indices32.push_back(baseIndex+i);
		meshData.Indices32.push_back(baseIndex+i+1);
	}

    return meshData;
}
 
void GeometryGenerator::Subdivide(MeshData& meshData)
{
	// Save a copy of the input geometry.
	MeshData inputCopy = meshData;


	meshData.Vertices.resize(0);
	meshData.Indices32.resize(0);

	//       v1
	//       *
	//      / \
	//     /   \
	//  m0*-----*m1
	//   / \   / \
	//  /   \ /   \
	// *-----*-----*
	// v0    m2     v2

	uint32 numTris = (uint32)inputCopy.Indices32.size()/3;
	for(uint32 i = 0; i < numTris; ++i)
	{
		Vertex v0 = inputCopy.Vertices[ inputCopy.Indices32[i*3+0] ];
		Vertex v1 = inputCopy.Vertices[ inputCopy.Indices32[i*3+1] ];
		Vertex v2 = inputCopy.Vertices[ inputCopy.Indices32[i*3+2] ];

		//
		// Generate the midpoints.
		//

        Vertex m0 = MidPoint(v0, v1);
        Vertex m1 = MidPoint(v1, v2);
        Vertex m2 = MidPoint(v0, v2);

		//
		// Add new geometry.
		//

		meshData.Vertices.push_back(v0); // 0
		meshData.Vertices.push_back(v1); // 1
		meshData.Vertices.push_back(v2); // 2
		meshData.Vertices.push_back(m0); // 3
		meshData.Vertices.push_back(m1); // 4
		meshData.Vertices.push_back(m2); // 5
 
		meshData.Indices32.push_back(i*6+0);
		meshData.Indices32.push_back(i*6+3);
		meshData.Indices32.push_back(i*6+5);

		meshData.Indices32.push_back(i*6+3);
		meshData.Indices32.push_back(i*6+4);
		meshData.Indices32.push_back(i*6+5);

		meshData.Indices32.push_back(i*6+5);
		meshData.Indices32.push_back(i*6+4);
		meshData.Indices32.push_back(i*6+2);

		meshData.Indices32.push_back(i*6+3);
		meshData.Indices32.push_back(i*6+1);
		meshData.Indices32.push_back(i*6+4);
	}
}

GeometryGenerator::Vertex GeometryGenerator::MidPoint(const Vertex& v0, const Vertex& v1)
{
    XMVECTOR p0 = XMLoadFloat3(&v0.Position);
    XMVECTOR p1 = XMLoadFloat3(&v1.Position);

    XMVECTOR n0 = XMLoadFloat3(&v0.Normal);
    XMVECTOR n1 = XMLoadFloat3(&v1.Normal);

    XMVECTOR tan0 = XMLoadFloat3(&v0.TangentU);
    XMVECTOR tan1 = XMLoadFloat3(&v1.TangentU);

    XMVECTOR tex0 = XMLoadFloat2(&v0.TexC);
    XMVECTOR tex1 = XMLoadFloat2(&v1.TexC);

    // Compute the midpoints of all the attributes.  Vectors need to be normalized
    // since linear interpolating can make them not unit length.  
    XMVECTOR pos = 0.5f*(p0 + p1);
    XMVECTOR normal = XMVector3Normalize(0.5f*(n0 + n1));
    XMVECTOR tangent = XMVector3Normalize(0.5f*(tan0+tan1));
    XMVECTOR tex = 0.5f*(tex0 + tex1);

    Vertex v;
    XMStoreFloat3(&v.Position, pos);
    XMStoreFloat3(&v.Normal, normal);
    XMStoreFloat3(&v.TangentU, tangent);
    XMStoreFloat2(&v.TexC, tex);

    return v;
}

GeometryGenerator::MeshData GeometryGenerator::CreateGeosphere(float radius, uint32 numSubdivisions)
{
    MeshData meshData;

	// Put a cap on the number of subdivisions.
    numSubdivisions = std::min<uint32>(numSubdivisions, 6u);

	// Approximate a sphere by tessellating an icosahedron.

	const float X = 0.525731f; 
	const float Z = 0.850651f;

	XMFLOAT3 pos[12] = 
	{
		XMFLOAT3(-X, 0.0f, Z),  XMFLOAT3(X, 0.0f, Z),  
		XMFLOAT3(-X, 0.0f, -Z), XMFLOAT3(X, 0.0f, -Z),    
		XMFLOAT3(0.0f, Z, X),   XMFLOAT3(0.0f, Z, -X), 
		XMFLOAT3(0.0f, -Z, X),  XMFLOAT3(0.0f, -Z, -X),    
		XMFLOAT3(Z, X, 0.0f),   XMFLOAT3(-Z, X, 0.0f), 
		XMFLOAT3(Z, -X, 0.0f),  XMFLOAT3(-Z, -X, 0.0f)
	};

    uint32 k[60] =
	{
		1,4,0,  4,9,0,  4,5,9,  8,5,4,  1,8,4,    
		1,10,8, 10,3,8, 8,3,5,  3,2,5,  3,7,2,    
		3,10,7, 10,6,7, 6,11,7, 6,0,11, 6,1,0, 
		10,1,6, 11,0,9, 2,11,9, 5,2,9,  11,2,7 
	};

    meshData.Vertices.resize(12);
    meshData.Indices32.assign(&k[0], &k[60]);

	for(uint32 i = 0; i < 12; ++i)
		meshData.Vertices[i].Position = pos[i];

	for(uint32 i = 0; i < numSubdivisions; ++i)
		Subdivide(meshData);

	// Project vertices onto sphere and scale.
	for(uint32 i = 0; i < meshData.Vertices.size(); ++i)
	{
		// Project onto unit sphere.
		XMVECTOR n = XMVector3Normalize(XMLoadFloat3(&meshData.Vertices[i].Position));

		// Project onto sphere.
		XMVECTOR p = radius*n;

		XMStoreFloat3(&meshData.Vertices[i].Position, p);
		XMStoreFloat3(&meshData.Vertices[i].Normal, n);

		// Derive texture coordinates from spherical coordinates.
        float theta = atan2f(meshData.Vertices[i].Position.z, meshData.Vertices[i].Position.x);

        // Put in [0, 2pi].
        if(theta < 0.0f)
            theta += XM_2PI;

		float phi = acosf(meshData.Vertices[i].Position.y / radius);

		meshData.Vertices[i].TexC.x = theta/XM_2PI;
		meshData.Vertices[i].TexC.y = phi/XM_PI;

		// Partial derivative of P with respect to theta
		meshData.Vertices[i].TangentU.x = -radius*sinf(phi)*sinf(theta);
		meshData.Vertices[i].TangentU.y = 0.0f;
		meshData.Vertices[i].TangentU.z = +radius*sinf(phi)*cosf(theta);

		XMVECTOR T = XMLoadFloat3(&meshData.Vertices[i].TangentU);
		XMStoreFloat3(&meshData.Vertices[i].TangentU, XMVector3Normalize(T));
	}

    return meshData;
}

GeometryGenerator::MeshData GeometryGenerator::CreateCylinder(float bottomRadius, float topRadius, float height, uint32 sliceCount, uint32 stackCount)
{
    MeshData meshData;

	//
	// Build Stacks.
	// 

	float stackHeight = height / stackCount;

	// Amount to increment radius as we move up each stack level from bottom to top.
	float radiusStep = (topRadius - bottomRadius) / stackCount;

	uint32 ringCount = stackCount+1;

	// Compute vertices for each stack ring starting at the bottom and moving up.
	for(uint32 i = 0; i < ringCount; ++i)
	{
		float y = -0.5f*height + i*stackHeight;
		float r = bottomRadius + i*radiusStep;

		// vertices of ring
		float dTheta = 2.0f*XM_PI/sliceCount;
		for(uint32 j = 0; j <= sliceCount; ++j)
		{
			Vertex vertex;

			float c = cosf(j*dTheta);
			float s = sinf(j*dTheta);

			vertex.Position = XMFLOAT3(r*c, y, r*s);

			vertex.TexC.x = (float)j/sliceCount;
			vertex.TexC.y = 1.0f - (float)i/stackCount;

			// Cylinder can be parameterized as follows, where we introduce v
			// parameter that goes in the same direction as the v tex-coord
			// so that the bitangent goes in the same direction as the v tex-coord.
			//   Let r0 be the bottom radius and let r1 be the top radius.
			//   y(v) = h - hv for v in [0,1].
			//   r(v) = r1 + (r0-r1)v
			//
			//   x(t, v) = r(v)*cos(t)
			//   y(t, v) = h - hv
			//   z(t, v) = r(v)*sin(t)
			// 
			//  dx/dt = -r(v)*sin(t)
			//  dy/dt = 0
			//  dz/dt = +r(v)*cos(t)
			//
			//  dx/dv = (r0-r1)*cos(t)
			//  dy/dv = -h
			//  dz/dv = (r0-r1)*sin(t)

			// This is unit length.
			vertex.TangentU = XMFLOAT3(-s, 0.0f, c);

			float dr = bottomRadius-topRadius;
			XMFLOAT3 bitangent(dr*c, -height, dr*s);

			XMVECTOR T = XMLoadFloat3(&vertex.TangentU);
			XMVECTOR B = XMLoadFloat3(&bitangent);
			XMVECTOR N = XMVector3Normalize(XMVector3Cross(T, B));
			XMStoreFloat3(&vertex.Normal, N);

			meshData.Vertices.push_back(vertex);
		}
	}

	// Add one because we duplicate the first and last vertex per ring
	// since the texture coordinates are different.
	uint32 ringVertexCount = sliceCount+1;

	// Compute indices for each stack.
	for(uint32 i = 0; i < stackCount; ++i)
	{
		for(uint32 j = 0; j < sliceCount; ++j)
		{
			meshData.Indices32.push_back(i*ringVertexCount + j);
			meshData.Indices32.push_back((i+1)*ringVertexCount + j);
			meshData.Indices32.push_back((i+1)*ringVertexCount + j+1);

			meshData.Indices32.push_back(i*ringVertexCount + j);
			meshData.Indices32.push_back((i+1)*ringVertexCount + j+1);
			meshData.Indices32.push_back(i*ringVertexCount + j+1);
		}
	}

	BuildCylinderTopCap(bottomRadius, topRadius, height, sliceCount, stackCount, meshData);
	BuildCylinderBottomCap(bottomRadius, topRadius, height, sliceCount, stackCount, meshData);

    return meshData;
}

void GeometryGenerator::BuildCylinderTopCap(float bottomRadius, float topRadius, float height,
											uint32 sliceCount, uint32 stackCount, MeshData& meshData)
{
	uint32 baseIndex = (uint32)meshData.Vertices.size();

	float y = 0.5f*height;
	float dTheta = 2.0f*XM_PI/sliceCount;

	// Duplicate cap ring vertices because the texture coordinates and normals differ.
	for(uint32 i = 0; i <= sliceCount; ++i)
	{
		float x = topRadius*cosf(i*dTheta);
		float z = topRadius*sinf(i*dTheta);

		// Scale down by the height to try and make top cap texture coord area
		// proportional to base.
		float u = x/height + 0.5f;
		float v = z/height + 0.5f;

		meshData.Vertices.push_back( Vertex(x, y, z, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, u, v) );
	}

	// Cap center vertex.
	meshData.Vertices.push_back( Vertex(0.0f, y, 0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.5f, 0.5f) );

	// Index of center vertex.
	uint32 centerIndex = (uint32)meshData.Vertices.size()-1;

	for(uint32 i = 0; i < sliceCount; ++i)
	{
		meshData.Indices32.push_back(centerIndex);
		meshData.Indices32.push_back(baseIndex + i+1);
		meshData.Indices32.push_back(baseIndex + i);
	}
}

void GeometryGenerator::BuildCylinderBottomCap(float bottomRadius, float topRadius, float height,
											   uint32 sliceCount, uint32 stackCount, MeshData& meshData)
{
	// 
	// Build bottom cap.
	//

	uint32 baseIndex = (uint32)meshData.Vertices.size();
	float y = -0.5f*height;

	// vertices of ring
	float dTheta = 2.0f*XM_PI/sliceCount;
	for(uint32 i = 0; i <= sliceCount; ++i)
	{
		float x = bottomRadius*cosf(i*dTheta);
		float z = bottomRadius*sinf(i*dTheta);

		// Scale down by the height to try and make top cap texture coord area
		// proportional to base.
		float u = x/height + 0.5f;
		float v = z/height + 0.5f;

		meshData.Vertices.push_back( Vertex(x, y, z, 0.0f, -1.0f, 0.0f, 1.0f, 0.0f, 0.0f, u, v) );
	}

	// Cap center vertex.
	meshData.Vertices.push_back( Vertex(0.0f, y, 0.0f, 0.0f, -1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.5f, 0.5f) );

	// Cache the index of center vertex.
	uint32 centerIndex = (uint32)meshData.Vertices.size()-1;

	for(uint32 i = 0; i < sliceCount; ++i)
	{
		meshData.Indices32.push_back(centerIndex);
		meshData.Indices32.push_back(baseIndex + i);
		meshData.Indices32.push_back(baseIndex + i+1);
	}
}

GeometryGenerator::MeshData GeometryGenerator::CreateGrid(float width, float depth, uint32 m, uint32 n)
{
    MeshData meshData;

	uint32 vertexCount = m*n;
	uint32 faceCount   = (m-1)*(n-1)*2;

	//
	// Create the vertices.
	//

	float halfWidth = 0.5f*width;
	float halfDepth = 0.5f*depth;

	float dx = width / (n-1);
	float dz = depth / (m-1);

	float du = 1.0f / (n-1);
	float dv = 1.0f / (m-1);

	meshData.Vertices.resize(vertexCount);
	for(uint32 i = 0; i < m; ++i)
	{
		float z = halfDepth - i*dz;
		for(uint32 j = 0; j < n; ++j)
		{
			float x = -halfWidth + j*dx;

			meshData.Vertices[i*n+j].Position = XMFLOAT3(x, 0.0f, z);
			meshData.Vertices[i*n+j].Normal   = XMFLOAT3(0.0f, 1.0f, 0.0f);
			meshData.Vertices[i*n+j].TangentU = XMFLOAT3(1.0f, 0.0f, 0.0f);

			// Stretch texture over grid.
			meshData.Vertices[i*n+j].TexC.x = j*du;
			meshData.Vertices[i*n+j].TexC.y = i*dv;
		}
	}
 
    //
	// Create the indices.
	//

	meshData.Indices32.resize(faceCount*3); // 3 indices per face

	// Iterate over each quad and compute indices.
	uint32 k = 0;
	for(uint32 i = 0; i < m-1; ++i)
	{
		for(uint32 j = 0; j < n-1; ++j)
		{
			meshData.Indices32[k]   = i*n+j;
			meshData.Indices32[k+1] = i*n+j+1;
			meshData.Indices32[k+2] = (i+1)*n+j;

			meshData.Indices32[k+3] = (i+1)*n+j;
			meshData.Indices32[k+4] = i*n+j+1;
			meshData.Indices32[k+5] = (i+1)*n+j+1;

			k += 6; // next quad
		}
	}

    return meshData;
}

GeometryGenerator::MeshData GeometryGenerator::CreateQuad(float x, float y, float w, float h, float depth)
{
    MeshData meshData;

	meshData.Vertices.resize(4);
	meshData.Indices32.resize(6);

	// Position coordinates specified in NDC space.
	meshData.Vertices[0] = Vertex(
        x, y - h, depth,
		0.0f, 0.0f, -1.0f,
		1.0f, 0.0f, 0.0f,
		0.0f, 1.0f);

	meshData.Vertices[1] = Vertex(
		x, y, depth,
		0.0f, 0.0f, -1.0f,
		1.0f, 0.0f, 0.0f,
		0.0f, 0.0f);

	meshData.Vertices[2] = Vertex(
		x+w, y, depth,
		0.0f, 0.0f, -1.0f,
		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f);

	meshData.Vertices[3] = Vertex(
		x+w, y-h, depth,
		0.0f, 0.0f, -1.0f,
		1.0f, 0.0f, 0.0f,
		1.0f, 1.0f);

	meshData.Indices32[0] = 0;
	meshData.Indices32[1] = 1;
	meshData.Indices32[2] = 2;

	meshData.Indices32[3] = 0;
	meshData.Indices32[4] = 2;
	meshData.Indices32[5] = 3;

    return meshData;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// Our solids start here!!
//////////////////////////////////////////////////////////////////////////////////////////////////

GeometryGenerator::MeshData GeometryGenerator::CreateStar()
{
	MeshData meshData;
	meshData.Vertices.resize(100);
	meshData.Indices32.resize(44*3);
	int i = 0;
	// FRONT FACE
	meshData.Vertices[i++] = Vertex(0.0f,     0.5f,   0.25f); // 0
	meshData.Vertices[i++] = Vertex(0.146f, 0.354f,   0.25f); // 1
	meshData.Vertices[i++] = Vertex(0.354f, 0.354f,   0.25f); // 2
	meshData.Vertices[i++] = Vertex(0.354f, 0.146f,   0.25f); // 3
	meshData.Vertices[i++] = Vertex(0.5f,   0.0f,     0.25f); // 4
	meshData.Vertices[i++] = Vertex(0.354f, -0.146f,  0.25f); // 5
	meshData.Vertices[i++] = Vertex(0.354f, -0.354f,  0.25f); // 6
	meshData.Vertices[i++] = Vertex(0.354f, -0.354f,  0.25f); // 7
	meshData.Vertices[i++] = Vertex(0.146f, -0.354f,  0.25f); // 8
	meshData.Vertices[i++] = Vertex(0.0f,     -0.5f,  0.25f); // 9
	meshData.Vertices[i++] = Vertex(-0.146f, -0.354f, 0.25f); // 10
	meshData.Vertices[i++] = Vertex(-0.354f, -0.354f, 0.25f); // 11
	meshData.Vertices[i++] = Vertex(-0.354f, -0.146f, 0.25f); // 12
	meshData.Vertices[i++] = Vertex(-0.5f,   0.0f,    0.25f); // 13
	meshData.Vertices[i++] = Vertex(-0.354f, 0.146f,  0.25f); // 14
	meshData.Vertices[i++] = Vertex(-0.354f, 0.354f,  0.25f); // 15
	meshData.Vertices[i++] = Vertex(-0.354f, 0.354f,  0.25f); // 16
	meshData.Vertices[i++] = Vertex(-0.146f, 0.354f,  0.25f); // 17
	// THE SIDES
	// NOTH EAST RECTANGLES
	meshData.Vertices[i++] = Vertex(0,     0.5f,     0.25f); // 18
	meshData.Vertices[i++] = Vertex(0,     0.5f,    -0.25f); // 19
	meshData.Vertices[i++] = Vertex(0.146f, 0.354f,  -0.25f); // 20
	meshData.Vertices[i++] = Vertex(0.146f, 0.354f,   0.25f); // 21
	meshData.Vertices[i++] = Vertex(0.146f, 0.354f,   0.25f); // 22
	meshData.Vertices[i++] = Vertex(0.146f, 0.354f,  -0.25f); // 23
	meshData.Vertices[i++] = Vertex(0.354f, 0.354f,  -0.25f); // 24
	meshData.Vertices[i++] = Vertex(0.354f, 0.354f,   0.25f); // 25
	meshData.Vertices[i++] = Vertex(0.354f, 0.354f,   0.25f); // 26
	meshData.Vertices[i++] = Vertex(0.354f, 0.354f,  -0.25f); // 27
	meshData.Vertices[i++] = Vertex(0.354f, 0.146f,   0.25f); // 28
	meshData.Vertices[i++] = Vertex(0.354f, 0.146f,  -0.25f); // 29
	meshData.Vertices[i++] = Vertex(0.354f, 0.146f,  -0.25f); // 30
	meshData.Vertices[i++] = Vertex(0.5f,   0,      -0.25f); // 31
	meshData.Vertices[i++] = Vertex(0.354f, 0.146f,   0.25f); // 32
	meshData.Vertices[i++] = Vertex(0.5f,   0,       0.25f); // 33

		// SOUTH EAST RECTANGLES
	meshData.Vertices[i++] = Vertex(0.5f,   0,       0.25f); // 34
	meshData.Vertices[i++] = Vertex(0.5f,   0,      -0.25f); // 35
	meshData.Vertices[i++] = Vertex(0.354f, -0.146f,  0.25f); // 36
	meshData.Vertices[i++] = Vertex(0.354f, -0.146f, -0.25f); // 37
	meshData.Vertices[i++] = Vertex(0.354f, -0.146f,  0.25f); // 38
	meshData.Vertices[i++] = Vertex(0.354f, -0.146f, -0.25f); // 39
	meshData.Vertices[i++] = Vertex(0.354f, -0.354f,  0.25f); // 40
	meshData.Vertices[i++] = Vertex(0.354f, -0.354f, -0.25f); // 41
	meshData.Vertices[i++] = Vertex(0.146f, -0.354f,  0.25f); // 42
	meshData.Vertices[i++] = Vertex(0.354f, -0.354f,  0.25f); // 43
	meshData.Vertices[i++] = Vertex(0.146f, -0.354f, -0.25f); // 44
	meshData.Vertices[i++] = Vertex(0.354f, -0.354f, -0.25f); // 45
	meshData.Vertices[i++] = Vertex(0.146f, -0.354f,  0.25f); // 46
	meshData.Vertices[i++] = Vertex(0.146f, -0.354f, -0.25f); // 47
	meshData.Vertices[i++] = Vertex(0,     -0.5f,    0.25f); // 48
	meshData.Vertices[i++] = Vertex(0,     -0.5f,   -0.25f); // 49

		// SOUTH WEST RECTANGLES
	meshData.Vertices[i++] = Vertex(-0.146f, -0.354f, 0.25f); // 50
	meshData.Vertices[i++] = Vertex(-0.146f, -0.354f,-0.25f); // 51
	meshData.Vertices[i++] = Vertex(0,     -0.5f,     0.25f); // 52
	meshData.Vertices[i++] = Vertex(0,     -0.5f,    -0.25f); // 53
	meshData.Vertices[i++] = Vertex(-0.354f, -0.354f, 0.25f); // 54
	meshData.Vertices[i++] = Vertex(-0.146f, -0.354f, 0.25f); // 55
	meshData.Vertices[i++] = Vertex(-0.354f, -0.354f,-0.25f); // 56
	meshData.Vertices[i++] = Vertex(-0.146f, -0.354f,-0.25f); // 57
	meshData.Vertices[i++] = Vertex(-0.354f, -0.146f,-0.25f); // 58
	meshData.Vertices[i++] = Vertex(-0.354f, -0.146f, 0.25f); // 59
	meshData.Vertices[i++] = Vertex(-0.354f, -0.354f,-0.25f); // 60
	meshData.Vertices[i++] = Vertex(-0.354f, -0.354f, 0.25f); // 61
	meshData.Vertices[i++] = Vertex(-0.5f,   0,       0.25f); // 62
	meshData.Vertices[i++] = Vertex(-0.354f, -0.146f, 0.25f); // 63
	meshData.Vertices[i++] = Vertex(-0.5f,   0,      -0.25f); // 64
	meshData.Vertices[i++] = Vertex(-0.354f, -0.146f,-0.25f); // 65

		// NORTH WEST RECTANGLES
	meshData.Vertices[i++] = Vertex(-0.354f, 0.146f, -0.25f); // 66
	meshData.Vertices[i++] = Vertex(-0.354f, 0.146f,  0.25f); // 67
	meshData.Vertices[i++] = Vertex(-0.5f,   0,      -0.25f); // 68
	meshData.Vertices[i++] = Vertex(-0.5f,   0,       0.25f); // 69
	meshData.Vertices[i++] = Vertex(-0.354f, 0.354f, -0.25f); // 70
	meshData.Vertices[i++] = Vertex(-0.354f, 0.354f,  0.25f); // 71
	meshData.Vertices[i++] = Vertex(-0.354f, 0.146f, -0.25f); // 72
	meshData.Vertices[i++] = Vertex(-0.354f, 0.146f,  0.25f); // 73
	meshData.Vertices[i++] = Vertex(-0.354f, 0.354f, -0.25f); // 74
	meshData.Vertices[i++] = Vertex(-0.146f, 0.354f, -0.25f); // 75
	meshData.Vertices[i++] = Vertex(-0.354f, 0.354f,  0.25f); // 76
	meshData.Vertices[i++] = Vertex(-0.146f, 0.354f,  0.25f); // 77
	meshData.Vertices[i++] = Vertex(-0.146f, 0.354f, -0.25f); // 78
	meshData.Vertices[i++] = Vertex(0,     0.5f,     -0.25f); // 79
	meshData.Vertices[i++] = Vertex(-0.146f, 0.354f,  0.25f); // 80
	meshData.Vertices[i++] = Vertex(0,     0.5f,      0.25f); // 81

		// BACK FACE
	meshData.Vertices[i++] = Vertex(0,     0.5f,     -0.25f); // 82
	meshData.Vertices[i++] = Vertex(0.146f, 0.354f,  -0.25f); // 83
	meshData.Vertices[i++] = Vertex(0.354f, 0.354f,  -0.25f); // 84
	meshData.Vertices[i++] = Vertex(0.354f, 0.146f,  -0.25f); // 85
	meshData.Vertices[i++] = Vertex(0.5f,   0,       -0.25f); // 86
	meshData.Vertices[i++] = Vertex(0.354f, -0.146f, -0.25f); // 87
	meshData.Vertices[i++] = Vertex(0.354f, -0.354f, -0.25f); // 88
	meshData.Vertices[i++] = Vertex(0.354f, -0.354f, -0.25f); // 89
	meshData.Vertices[i++] = Vertex(0.146f, -0.354f, -0.25f); // 90
	meshData.Vertices[i++] = Vertex(0,     -0.5f,    -0.25f); // 91
	meshData.Vertices[i++] = Vertex(-0.146f, -0.354f,-0.25f); // 92
	meshData.Vertices[i++] = Vertex(-0.354f, -0.354f,-0.25f); // 93
	meshData.Vertices[i++] = Vertex(-0.354f, -0.146f,-0.25f); // 94
	meshData.Vertices[i++] = Vertex(-0.5f,   0,      -0.25f); // 95
	meshData.Vertices[i++] = Vertex(-0.354f, 0.146f, -0.25f); // 96
	meshData.Vertices[i++] = Vertex(-0.354f, 0.354f, -0.25f); // 97
	meshData.Vertices[i++] = Vertex(-0.354f, 0.354f, -0.25f); // 98
	meshData.Vertices[i++] = Vertex(-0.146f, 0.354f, -0.25f); // 99

	i = 0; // Restart numbering
	meshData.Indices32[i++] = 6;  meshData.Indices32[i++] = 2; meshData.Indices32[i++] = 16; 
	meshData.Indices32[i++] = 11; meshData.Indices32[i++] = 7; meshData.Indices32[i++] = 15;
	meshData.Indices32[i++] = 17; meshData.Indices32[i++] = 1; meshData.Indices32[i++] = 0;
	meshData.Indices32[i++] = 5;  meshData.Indices32[i++] = 4; meshData.Indices32[i++] = 3;
	meshData.Indices32[i++] = 10; meshData.Indices32[i++] = 9; meshData.Indices32[i++] = 8;
	meshData.Indices32[i++] = 14; meshData.Indices32[i++] = 13; meshData.Indices32[i++] = 12;
	meshData.Indices32[i++] = 21; meshData.Indices32[i++] = 20; meshData.Indices32[i++] = 19;
	meshData.Indices32[i++] = 19; meshData.Indices32[i++] = 18; meshData.Indices32[i++] = 21;
	meshData.Indices32[i++] = 25; meshData.Indices32[i++] = 24; meshData.Indices32[i++] = 23;
	meshData.Indices32[i++] = 23; meshData.Indices32[i++] = 22; meshData.Indices32[i++] = 25;
	meshData.Indices32[i++] = 29; meshData.Indices32[i++] = 27; meshData.Indices32[i++] = 26;
	meshData.Indices32[i++] = 26; meshData.Indices32[i++] = 28; meshData.Indices32[i++] = 29;
	meshData.Indices32[i++] = 33; meshData.Indices32[i++] = 31; meshData.Indices32[i++] = 30;
	meshData.Indices32[i++] = 30; meshData.Indices32[i++] = 32; meshData.Indices32[i++] = 33;
	meshData.Indices32[i++] = 37; meshData.Indices32[i++] = 35; meshData.Indices32[i++] = 34;
	meshData.Indices32[i++] = 34; meshData.Indices32[i++] = 36; meshData.Indices32[i++] = 37;
	meshData.Indices32[i++] = 41; meshData.Indices32[i++] = 39; meshData.Indices32[i++] = 38;
	meshData.Indices32[i++] = 38; meshData.Indices32[i++] = 40; meshData.Indices32[i++] = 41;
	meshData.Indices32[i++] = 45; meshData.Indices32[i++] = 43; meshData.Indices32[i++] = 42;
	meshData.Indices32[i++] = 42; meshData.Indices32[i++] = 44; meshData.Indices32[i++] = 45;
	meshData.Indices32[i++] = 47; meshData.Indices32[i++] = 46; meshData.Indices32[i++] = 48;
	meshData.Indices32[i++] = 48; meshData.Indices32[i++] = 49; meshData.Indices32[i++] = 47;
	meshData.Indices32[i++] = 53; meshData.Indices32[i++] = 52; meshData.Indices32[i++] = 50;
	meshData.Indices32[i++] = 50; meshData.Indices32[i++] = 51; meshData.Indices32[i++] = 53;
	meshData.Indices32[i++] = 57; meshData.Indices32[i++] = 55; meshData.Indices32[i++] = 54;
	meshData.Indices32[i++] = 54; meshData.Indices32[i++] = 56; meshData.Indices32[i++] = 57;
	meshData.Indices32[i++] = 61; meshData.Indices32[i++] = 59; meshData.Indices32[i++] = 58;
	meshData.Indices32[i++] = 58; meshData.Indices32[i++] = 60; meshData.Indices32[i++] = 61;
	meshData.Indices32[i++] = 65; meshData.Indices32[i++] = 63; meshData.Indices32[i++] = 62;
	meshData.Indices32[i++] = 62; meshData.Indices32[i++] = 64; meshData.Indices32[i++] = 65;
	meshData.Indices32[i++] = 67; meshData.Indices32[i++] = 66; meshData.Indices32[i++] = 68;
	meshData.Indices32[i++] = 68; meshData.Indices32[i++] = 69; meshData.Indices32[i++] = 67;
	meshData.Indices32[i++] = 73; meshData.Indices32[i++] = 71; meshData.Indices32[i++] = 70;
	meshData.Indices32[i++] = 70; meshData.Indices32[i++] = 72; meshData.Indices32[i++] = 73;
	meshData.Indices32[i++] = 77; meshData.Indices32[i++] = 75; meshData.Indices32[i++] = 74;
	meshData.Indices32[i++] = 74; meshData.Indices32[i++] = 76; meshData.Indices32[i++] = 77;
	meshData.Indices32[i++] = 81; meshData.Indices32[i++] = 79; meshData.Indices32[i++] = 78;
	meshData.Indices32[i++] = 78; meshData.Indices32[i++] = 80; meshData.Indices32[i++] = 81;
	meshData.Indices32[i++] = 98; meshData.Indices32[i++] = 84; meshData.Indices32[i++] = 88;
	meshData.Indices32[i++] = 93; meshData.Indices32[i++] = 97; meshData.Indices32[i++] = 89;
	meshData.Indices32[i++] = 99; meshData.Indices32[i++] = 82; meshData.Indices32[i++] = 83;
	meshData.Indices32[i++] = 85; meshData.Indices32[i++] = 86; meshData.Indices32[i++] = 87;
	meshData.Indices32[i++] = 91; meshData.Indices32[i++] = 92; meshData.Indices32[i++] = 90;
	meshData.Indices32[i++] = 95; meshData.Indices32[i++] = 96; meshData.Indices32[i++] = 94;

	return meshData;

}


GeometryGenerator::MeshData GeometryGenerator::CreateWedge()
{
	MeshData meshData;

	meshData.Vertices.resize(18);
	meshData.Indices32.resize(8 * 3);
	int i = 0;
	meshData.Vertices[i++] = Vertex(0,  1, -1); // 0
	meshData.Vertices[i++] = Vertex(0,  1,  0); // 1
	meshData.Vertices[i++] = Vertex(0,  0, -1); // 2
	meshData.Vertices[i++] = Vertex(0,  0,  0); // 3
	meshData.Vertices[i++] = Vertex(0,  1,  0); // 4
	meshData.Vertices[i++] = Vertex(0,  0,  0); // 5
	meshData.Vertices[i++] = Vertex(1,  0,  0); // 6
	meshData.Vertices[i++] = Vertex(0,  0,  0); // 7
	meshData.Vertices[i++] = Vertex(1,  0,  0); // 8
	meshData.Vertices[i++] = Vertex(0,  0, -1); // 9
	meshData.Vertices[i++] = Vertex(1,  0, -1); // 10
	meshData.Vertices[i++] = Vertex(0,  1,  0); // 11
	meshData.Vertices[i++] = Vertex(0,  1, -1); // 12
	meshData.Vertices[i++] = Vertex(1,  0,  0); // 13
	meshData.Vertices[i++] = Vertex(1,  0, -1); // 14
	meshData.Vertices[i++] = Vertex(0,  1, -1); // 15
	meshData.Vertices[i++] = Vertex(0,  0, -1); // 16
	meshData.Vertices[i++] = Vertex(1,  0, -1); // 17

	// Centering X, Y and Z to account for transformations
	
	for (int j = 0; j < 18; j++) {
		meshData.Vertices[j].Position.x -= 0.5f;
		meshData.Vertices[j].Position.y -= 0.5f;
		meshData.Vertices[j].Position.z += 0.5f;
	}
	i = 0;
	meshData.Indices32[i++] = 3;  meshData.Indices32[i++] = 1;  meshData.Indices32[i++] = 0;
	meshData.Indices32[i++] = 0;  meshData.Indices32[i++] = 2;  meshData.Indices32[i++] = 3;
	meshData.Indices32[i++] = 5;  meshData.Indices32[i++] = 6;  meshData.Indices32[i++] = 4;
	meshData.Indices32[i++] = 10; meshData.Indices32[i++] = 8;  meshData.Indices32[i++] = 7;
	meshData.Indices32[i++] = 7;  meshData.Indices32[i++] = 9;  meshData.Indices32[i++] = 10;
	meshData.Indices32[i++] = 13; meshData.Indices32[i++] = 14; meshData.Indices32[i++] = 12;
	meshData.Indices32[i++] = 12; meshData.Indices32[i++] = 11; meshData.Indices32[i++] = 13;
	meshData.Indices32[i++] = 17; meshData.Indices32[i++] = 16; meshData.Indices32[i++] = 15;

	return meshData;
}


GeometryGenerator::MeshData  GeometryGenerator::CreateTriPrism()
{
	MeshData meshData;

	meshData.Vertices.resize(12);
	meshData.Indices32.resize(4 * 3);
	int i = 0;
	meshData.Vertices[i++] = Vertex(0,     0,        0     );// 0
	meshData.Vertices[i++] = Vertex(0.5f,  0.866f,  -0.28f );// 1
	meshData.Vertices[i++] = Vertex(1,     0,        0     );// 2
	meshData.Vertices[i++] = Vertex(0,     0,        0     );// 3
	meshData.Vertices[i++] = Vertex(1,     0,        0     );// 4
	meshData.Vertices[i++] = Vertex(0.5f,  0,       -0.866f);// 5
	meshData.Vertices[i++] = Vertex(0.5f,  0.866f,  -0.289f);// 6
	meshData.Vertices[i++] = Vertex(0.5f,  0,       -0.866f);// 7
	meshData.Vertices[i++] = Vertex(1,     0,        0     );// 8
	meshData.Vertices[i++] = Vertex(0.5f,  0,       -0.866f);// 9
	meshData.Vertices[i++] = Vertex(0.5f,  0.866f,  -0.289f);// 10
	meshData.Vertices[i++] = Vertex(0,     0,        0     );// 11

	// Centering X and Y
	
	for (int j = 0; j < 12; j++) {
		meshData.Vertices[j].Position.x -= 0.5f;
		meshData.Vertices[j].Position.y -= 0.5f;
		meshData.Vertices[j].Position.z += 0.5f;
	}
	i = 0;
	meshData.Indices32[i++] = 2; meshData.Indices32[i++] = 1; meshData.Indices32[i++] = 0;
	meshData.Indices32[i++] = 5; meshData.Indices32[i++] = 4; meshData.Indices32[i++] = 3;
	meshData.Indices32[i++] = 8; meshData.Indices32[i++] = 7; meshData.Indices32[i++] = 6;
	meshData.Indices32[i++] = 11; meshData.Indices32[i++] = 10; meshData.Indices32[i++] = 9;

	return meshData;
}

GeometryGenerator::MeshData GeometryGenerator::CreateCube() {

	MeshData meshData;

	// CUBE VERTICES VBO -----------------------------------------------------------
	float halfCubeSide = 0.5f;
	meshData.Vertices.resize(24);
	meshData.Indices32.resize(12 * 3);
	int i = 0;
	/*
	 * This is the first set - from 0 to 7
	 */
	meshData.Vertices[i++] = Vertex(-halfCubeSide, -halfCubeSide,  halfCubeSide);	// 0,1                            
	meshData.Vertices[i++] = Vertex( halfCubeSide, -halfCubeSide,  halfCubeSide);	// 1,1
	meshData.Vertices[i++] = Vertex( halfCubeSide,  halfCubeSide,  halfCubeSide);	// 1,0
	meshData.Vertices[i++] = Vertex(-halfCubeSide,  halfCubeSide,  halfCubeSide);	// 0,0
	meshData.Vertices[i++] = Vertex(-halfCubeSide, -halfCubeSide, -halfCubeSide);	// 0,1
	meshData.Vertices[i++] = Vertex( halfCubeSide, -halfCubeSide, -halfCubeSide);	// 1,1
	meshData.Vertices[i++] = Vertex( halfCubeSide,  halfCubeSide, -halfCubeSide);	// 1,0
	meshData.Vertices[i++] = Vertex(-halfCubeSide,  halfCubeSide, -halfCubeSide);	// 0,0
	/*
 	 * This is the second set - from 8 to 15
	 */
	meshData.Vertices[i++] = Vertex(-halfCubeSide, -halfCubeSide,  halfCubeSide);	// 0,1
	meshData.Vertices[i++] = Vertex( halfCubeSide, -halfCubeSide,  halfCubeSide);	// 1,1
	meshData.Vertices[i++] = Vertex( halfCubeSide,  halfCubeSide,  halfCubeSide);	// 1,0
	meshData.Vertices[i++] = Vertex(-halfCubeSide,  halfCubeSide,  halfCubeSide);	// 0,0
	meshData.Vertices[i++] = Vertex(-halfCubeSide, -halfCubeSide, -halfCubeSide);	// 0,1
	meshData.Vertices[i++] = Vertex( halfCubeSide, -halfCubeSide, -halfCubeSide);	// 1,1
	meshData.Vertices[i++] = Vertex( halfCubeSide,  halfCubeSide, -halfCubeSide);	// 1,0
	meshData.Vertices[i++] = Vertex(-halfCubeSide,  halfCubeSide, -halfCubeSide);	// 0,0
	/*
	 * This is the third set - from 16 to 23
	 */
	meshData.Vertices[i++] = Vertex(-halfCubeSide, -halfCubeSide,  halfCubeSide);	// 0,1
	meshData.Vertices[i++] = Vertex( halfCubeSide, -halfCubeSide,  halfCubeSide);	// 1,1
	meshData.Vertices[i++] = Vertex( halfCubeSide,  halfCubeSide,  halfCubeSide);	// 1,0
	meshData.Vertices[i++] = Vertex(-halfCubeSide,  halfCubeSide,  halfCubeSide);	// 0,0
	meshData.Vertices[i++] = Vertex(-halfCubeSide, -halfCubeSide, -halfCubeSide);	// 0,1
	meshData.Vertices[i++] = Vertex( halfCubeSide, -halfCubeSide, -halfCubeSide);	// 1,1
	meshData.Vertices[i++] = Vertex( halfCubeSide,  halfCubeSide, -halfCubeSide);	// 1,0
	meshData.Vertices[i++] = Vertex(-halfCubeSide,  halfCubeSide, -halfCubeSide);	// 0,0

	i = 0;
	// front
	meshData.Indices32[i++] = 0;  meshData.Indices32[i++] = 1;   meshData.Indices32[i++] = 2;
	meshData.Indices32[i++] = 2;  meshData.Indices32[i++] = 3;   meshData.Indices32[i++] = 0;
	// back
	meshData.Indices32[i++] = 4;  meshData.Indices32[i++] = 7;   meshData.Indices32[i++] = 6;
	meshData.Indices32[i++] = 6;  meshData.Indices32[i++] = 5;   meshData.Indices32[i++] = 4;
	// top
	meshData.Indices32[i++] = 11; meshData.Indices32[i++] = 10;  meshData.Indices32[i++] = 14;
	meshData.Indices32[i++] = 14; meshData.Indices32[i++] = 15;  meshData.Indices32[i++] = 11;
	// bottom
	meshData.Indices32[i++] = 8;  meshData.Indices32[i++] = 12;  meshData.Indices32[i++] = 13;
	meshData.Indices32[i++] = 13; meshData.Indices32[i++] = 9;   meshData.Indices32[i++] = 8;
	// right
	meshData.Indices32[i++] = 17; meshData.Indices32[i++] = 21;  meshData.Indices32[i++] = 22;
	meshData.Indices32[i++] = 22; meshData.Indices32[i++] = 18;  meshData.Indices32[i++] = 17;
	// left
	meshData.Indices32[i++] = 19; meshData.Indices32[i++] =  23; meshData.Indices32[i++] =  20;
	meshData.Indices32[i++] = 20; meshData.Indices32[i++] = 16; meshData.Indices32[i++] = 19;

	return meshData;
}

GeometryGenerator::MeshData GeometryGenerator::CreatePyramid() {
	MeshData meshData;
	meshData.Vertices.resize(16);
	meshData.Indices32.resize(6 * 3);
	int i = 0;

	// PYRAMID VERTICES VBO -----------------------------------------------------------
	float halfPyramidSide = 0.5f;
	float pyramidHeight = 1.0f;

	meshData.Vertices[i++] = Vertex(-halfPyramidSide, 0.0f,  halfPyramidSide); // 00-FRONT_1
	meshData.Vertices[i++] = Vertex( halfPyramidSide, 0.0f,  halfPyramidSide); // 01-FRONT_2
	meshData.Vertices[i++] = Vertex( halfPyramidSide, 0.0f, -halfPyramidSide); // 02-BACK_1
	meshData.Vertices[i++] = Vertex(-halfPyramidSide, 0.0f, -halfPyramidSide); // 03-BACK_2
	meshData.Vertices[i++] = Vertex( 0.0f,            pyramidHeight, 0.0f	); // 04-TOP_1
	meshData.Vertices[i++] = Vertex(-halfPyramidSide, 0.0f,  halfPyramidSide); // 05-LEFT_2
	meshData.Vertices[i++] = Vertex( halfPyramidSide, 0.0f,  halfPyramidSide); // 06-RIGHT_1
	meshData.Vertices[i++] = Vertex( halfPyramidSide, 0.0f, -halfPyramidSide); // 07-RIGHT_2
	meshData.Vertices[i++] = Vertex(-halfPyramidSide, 0.0f, -halfPyramidSide); // 08-LEFT_1
	meshData.Vertices[i++] = Vertex( 0.0f,            pyramidHeight, 0.0f	); // 09-TOP_2
	meshData.Vertices[i++] = Vertex(-halfPyramidSide, 0.0f,  halfPyramidSide); // 10-BOTTOM_1
	meshData.Vertices[i++] = Vertex( halfPyramidSide, 0.0f,  halfPyramidSide); // 11-BOTTOM_2
	meshData.Vertices[i++] = Vertex( halfPyramidSide, 0.0f, -halfPyramidSide); // 12-BOTTOM_3
	meshData.Vertices[i++] = Vertex(-halfPyramidSide, 0.0f, -halfPyramidSide); // 13-BOTTOM_4
	meshData.Vertices[i++] = Vertex( 0.0f,            pyramidHeight, 0.0f	); // 04-TOP_3
	meshData.Vertices[i++] = Vertex( 0.0f,            pyramidHeight, 0.0f	); // 04-TOP_4

	/*
	 Centering the pyramid
	*/
	for (int j = 0; j < 16; j++) {
		meshData.Vertices[j].Position.y -= pyramidHeight / 2;
	}
	i = 0;

		// front
	meshData.Indices32[i++] = 0;  meshData.Indices32[i++] = 1;  meshData.Indices32[i++] = 4;
	// back
	meshData.Indices32[i++] = 2;  meshData.Indices32[i++] = 3;  meshData.Indices32[i++] = 14;
	// left
	meshData.Indices32[i++] = 8;  meshData.Indices32[i++] = 5;  meshData.Indices32[i++] = 9;
	// right
	meshData.Indices32[i++] = 6;  meshData.Indices32[i++] = 7;  meshData.Indices32[i++] = 15;
	// bottom
	meshData.Indices32[i++] = 10; meshData.Indices32[i++] = 13; meshData.Indices32[i++] = 12;
	meshData.Indices32[i++] = 12; meshData.Indices32[i++] = 11; meshData.Indices32[i++] = 10;

	return meshData;

}

GeometryGenerator::MeshData GeometryGenerator::CreateHexagon() {

	MeshData meshData;
	meshData.Vertices.resize(24);
	meshData.Indices32.resize(20 * 3);
	int i = 0;

	// HEXAGON VERTICES VBO -----------------------------------------------------------
	float halfSideSize = 0.5f;
	float sdDist = halfSideSize * sqrt(3) / 2; // side * sin(60)

		// TOP_LID
	meshData.Vertices[i++] = Vertex(-halfSideSize,			halfSideSize,  halfSideSize);
	meshData.Vertices[i++] = Vertex( halfSideSize,			halfSideSize,  halfSideSize);
	meshData.Vertices[i++] = Vertex( halfSideSize + sdDist,	halfSideSize,  0		   );
	meshData.Vertices[i++] = Vertex( halfSideSize,			halfSideSize, -halfSideSize);
	meshData.Vertices[i++] = Vertex(-halfSideSize,			halfSideSize, -halfSideSize);
	meshData.Vertices[i++] = Vertex(-halfSideSize - sdDist,	halfSideSize,  0		   );

		// BOTTOM_LID
	meshData.Vertices[i++] = Vertex(-halfSideSize,		   -halfSideSize,  halfSideSize);
	meshData.Vertices[i++] = Vertex( halfSideSize,		   -halfSideSize,  halfSideSize);
	meshData.Vertices[i++] = Vertex( halfSideSize + sdDist,-halfSideSize,  0		   );
	meshData.Vertices[i++] = Vertex( halfSideSize,		   -halfSideSize, -halfSideSize);
	meshData.Vertices[i++] = Vertex(-halfSideSize,		   -halfSideSize, -halfSideSize);
	meshData.Vertices[i++] = Vertex(-halfSideSize - sdDist,-halfSideSize,  0		   );

	   // SIDES
	meshData.Vertices[i++] = Vertex(-halfSideSize,			halfSideSize,  halfSideSize);
	meshData.Vertices[i++] = Vertex(-halfSideSize,		   -halfSideSize,  halfSideSize);
	meshData.Vertices[i++] = Vertex( halfSideSize,			halfSideSize,  halfSideSize);
	meshData.Vertices[i++] = Vertex( halfSideSize,		   -halfSideSize,  halfSideSize);
	meshData.Vertices[i++] = Vertex( halfSideSize + sdDist,	halfSideSize,  0		   );
	meshData.Vertices[i++] = Vertex( halfSideSize + sdDist,-halfSideSize,  0		   );
	meshData.Vertices[i++] = Vertex( halfSideSize,			halfSideSize, -halfSideSize);
	meshData.Vertices[i++] = Vertex( halfSideSize,		   -halfSideSize, -halfSideSize);
	meshData.Vertices[i++] = Vertex(-halfSideSize,			halfSideSize, -halfSideSize);
	meshData.Vertices[i++] = Vertex(-halfSideSize,		   -halfSideSize, -halfSideSize);
	meshData.Vertices[i++] = Vertex(-halfSideSize - sdDist,	halfSideSize,  0		   );
	meshData.Vertices[i++] = Vertex(-halfSideSize - sdDist,-halfSideSize,  0		   );

	i = 0;
		// top
	meshData.Indices32[i++] = 0;  meshData.Indices32[i++] = 1;  meshData.Indices32[i++] = 3;
	meshData.Indices32[i++] = 1;  meshData.Indices32[i++] = 2;  meshData.Indices32[i++] = 3;
	meshData.Indices32[i++] = 0;  meshData.Indices32[i++] = 3;  meshData.Indices32[i++] = 5;
	meshData.Indices32[i++] = 3;  meshData.Indices32[i++] = 4;  meshData.Indices32[i++] = 5;
	//bottom					  							   
	meshData.Indices32[i++] = 0;  meshData.Indices32[i++] = 3;  meshData.Indices32[i++] = 1;
	meshData.Indices32[i++] = 1;  meshData.Indices32[i++] = 3;  meshData.Indices32[i++] = 2;
	meshData.Indices32[i++] = 5;  meshData.Indices32[i++] = 3;  meshData.Indices32[i++] = 0;
	meshData.Indices32[i++] = 5;  meshData.Indices32[i++] = 4;  meshData.Indices32[i++] = 3;
	// sides					  							   
	meshData.Indices32[i++] = 0;  meshData.Indices32[i++] = 1;  meshData.Indices32[i++] = 3;
	meshData.Indices32[i++] = 3;  meshData.Indices32[i++] = 2;  meshData.Indices32[i++] = 0;
	meshData.Indices32[i++] = 2;  meshData.Indices32[i++] = 3;  meshData.Indices32[i++] = 5;
	meshData.Indices32[i++] = 5;  meshData.Indices32[i++] = 4;  meshData.Indices32[i++] = 2;
	meshData.Indices32[i++] = 4;  meshData.Indices32[i++] = 5;  meshData.Indices32[i++] = 7;
	meshData.Indices32[i++] = 7;  meshData.Indices32[i++] = 6;  meshData.Indices32[i++] = 4;
	meshData.Indices32[i++] = 6;  meshData.Indices32[i++] = 7;  meshData.Indices32[i++] = 9;
	meshData.Indices32[i++] = 9;  meshData.Indices32[i++] = 8;  meshData.Indices32[i++] = 6;
	meshData.Indices32[i++] = 8;  meshData.Indices32[i++] = 9;  meshData.Indices32[i++] = 11;
	meshData.Indices32[i++] = 11; meshData.Indices32[i++] = 10; meshData.Indices32[i++] = 8;
	meshData.Indices32[i++] = 10; meshData.Indices32[i++] = 11; meshData.Indices32[i++] = 1;
	meshData.Indices32[i++] = 1;  meshData.Indices32[i++] =  0; meshData.Indices32[i++] =  10;

	// Shifting stuff
	for (int j = 4 * 3; j < 8 * 3; j++) {
		meshData.Indices32[j] += 6;
	}

	for (int j = 8 * 3; j < 20*3; j++) {
		meshData.Indices32[j] += 12;
	}


	return meshData;
}

GeometryGenerator::MeshData GeometryGenerator::CreateSimsIndicator() {

	MeshData meshData;
	meshData.Vertices.resize(36);
	meshData.Indices32.resize(20 * 3);
	int i = 0;

	// SIMS INDICATOR VERTICES VBO -----------------------------------------------------------
	float halfSideSize = 0.5f;
	float sdDist = halfSideSize * sqrt(3) / 2; // side * sin(60)

	meshData.Vertices[i++] = Vertex( 0,                     halfSideSize, 0           );
	meshData.Vertices[i++] = Vertex(-halfSideSize,			0,            halfSideSize);
	meshData.Vertices[i++] = Vertex( halfSideSize,			0,            halfSideSize);
	meshData.Vertices[i++] = Vertex( 0,                     halfSideSize, 0			  );
	meshData.Vertices[i++] = Vertex( halfSideSize,			0,            halfSideSize);
	meshData.Vertices[i++] = Vertex( halfSideSize + sdDist,	0,            0			  );
	meshData.Vertices[i++] = Vertex( 0,                     halfSideSize, 0			  );
	meshData.Vertices[i++] = Vertex( halfSideSize + sdDist,	0,            0			  );
	meshData.Vertices[i++] = Vertex( halfSideSize,			0,           -halfSideSize);
	meshData.Vertices[i++] = Vertex( 0,                     halfSideSize, 0			  );
	meshData.Vertices[i++] = Vertex( halfSideSize,			0,           -halfSideSize);
	meshData.Vertices[i++] = Vertex(-halfSideSize,			0,           -halfSideSize);
	meshData.Vertices[i++] = Vertex( 0,                     halfSideSize, 0			  );
	meshData.Vertices[i++] = Vertex(-halfSideSize,			0,           -halfSideSize);
	meshData.Vertices[i++] = Vertex(-halfSideSize - sdDist,	0,            0			  );
	meshData.Vertices[i++] = Vertex( 0,                     halfSideSize, 0			  );
	meshData.Vertices[i++] = Vertex(-halfSideSize - sdDist,	0,            0			  );
	meshData.Vertices[i++] = Vertex(-halfSideSize,			0,            halfSideSize);
	meshData.Vertices[i++] = Vertex( 0,                    -halfSideSize, 0			  );
	meshData.Vertices[i++] = Vertex( halfSideSize,			0,            halfSideSize);
	meshData.Vertices[i++] = Vertex(-halfSideSize,			0,            halfSideSize);
	meshData.Vertices[i++] = Vertex( 0,                    -halfSideSize, 0			  );
	meshData.Vertices[i++] = Vertex( halfSideSize + sdDist,	0,            0			  );
	meshData.Vertices[i++] = Vertex( halfSideSize,			0,            halfSideSize);
	meshData.Vertices[i++] = Vertex( 0,                    -halfSideSize, 0			  );
	meshData.Vertices[i++] = Vertex( halfSideSize,			0,           -halfSideSize);
	meshData.Vertices[i++] = Vertex( halfSideSize + sdDist,	0,            0			  );
	meshData.Vertices[i++] = Vertex( 0,                    -halfSideSize, 0			  );
	meshData.Vertices[i++] = Vertex(-halfSideSize,			0,           -halfSideSize);
	meshData.Vertices[i++] = Vertex( halfSideSize,			0,           -halfSideSize);
	meshData.Vertices[i++] = Vertex( 0,                    -halfSideSize, 0			  );
	meshData.Vertices[i++] = Vertex(-halfSideSize - sdDist,	0,            0			  );
	meshData.Vertices[i++] = Vertex(-halfSideSize,			0,           -halfSideSize);
	meshData.Vertices[i++] = Vertex( 0,                    -halfSideSize, 0			  );
	meshData.Vertices[i++] = Vertex(-halfSideSize,			0,            halfSideSize);
	meshData.Vertices[i++] = Vertex(-halfSideSize - sdDist,	0,            0			  );

	i = 0;
	meshData.Indices32[i++] = 0;  meshData.Indices32[i++] = 1;   meshData.Indices32[i++] = 2;
	meshData.Indices32[i++] = 3;  meshData.Indices32[i++] = 4;   meshData.Indices32[i++] = 5;
	meshData.Indices32[i++] = 6;  meshData.Indices32[i++] = 7;   meshData.Indices32[i++] = 8;
	meshData.Indices32[i++] = 9;  meshData.Indices32[i++] = 10;  meshData.Indices32[i++] = 11;
	meshData.Indices32[i++] = 12; meshData.Indices32[i++] = 13;  meshData.Indices32[i++] = 14;
	meshData.Indices32[i++] = 15; meshData.Indices32[i++] = 16;  meshData.Indices32[i++] = 17;
	meshData.Indices32[i++] = 18; meshData.Indices32[i++] = 19;  meshData.Indices32[i++] = 20;
	meshData.Indices32[i++] = 21; meshData.Indices32[i++] = 22;  meshData.Indices32[i++] = 23;
	meshData.Indices32[i++] = 24; meshData.Indices32[i++] = 25;  meshData.Indices32[i++] = 26;
	meshData.Indices32[i++] = 27; meshData.Indices32[i++] = 28;  meshData.Indices32[i++] = 29;
	meshData.Indices32[i++] = 30; meshData.Indices32[i++] =  31; meshData.Indices32[i++] = 32;
	meshData.Indices32[i++] = 33; meshData.Indices32[i++] = 34; meshData.Indices32[i++] = 35;

	return meshData;
}

GeometryGenerator::MeshData GeometryGenerator::CreateMyCylinder() {

	MeshData meshData;
	int i = 0;

	float radius = 0.5f;
	int degrees = 2;
	float halfHeight = 0.5f;

	int baseTriangles = 360 / degrees;
	int triangles = baseTriangles * 2 * 2;
	int vertices = triangles * 3;
	meshData.Vertices.resize(vertices);
	meshData.Indices32.resize(vertices);

	for (int currentAngle = 0; currentAngle < 360; currentAngle += degrees) {

		float x1 = radius * cos(currentAngle * PI / 180);
		float x2 = radius * cos((currentAngle + degrees) * PI / 180);
		float z1 = radius * sin(currentAngle * PI / 180);
		float z2 = radius * sin((currentAngle + degrees) * PI / 180);

		/* TOP Origin */
		meshData.Vertices[i++] = Vertex(0.0f, halfHeight, 0.0f);
		/* A */
		meshData.Vertices[i++] = Vertex(x1, halfHeight, z1);
		/* B */
		meshData.Vertices[i++] = Vertex(x2,halfHeight,z2);

		/* DOWN Origin */
		meshData.Vertices[i++] = Vertex(0.0f,-halfHeight, 0.0f);
		/* A */
		meshData.Vertices[i++] = Vertex(x2,-halfHeight,z2);
		/* B */
		meshData.Vertices[i++] = Vertex(x1,-halfHeight, z1);

		/* SIDE 1 */
		/* A */
		meshData.Vertices[i++] = Vertex(x2,halfHeight,z2);
		/* B */
		meshData.Vertices[i++] = Vertex(x1,halfHeight, z1);
		/* C */
		meshData.Vertices[i++] = Vertex(x1,-halfHeight, z1);

		/* SIDE 2 */
		/* C */
		meshData.Vertices[i++] = Vertex(x1,-halfHeight,z1); 
		/* D */										 
		meshData.Vertices[i++] = Vertex(x2,-halfHeight,z2); 
		/* A */										 
		meshData.Vertices[i++] = Vertex(x2, halfHeight,z2); 
	}

	for (int j = 0; j < vertices; j++) {
		meshData.Indices32[j] = j;
	}

	return meshData;
}

GeometryGenerator::MeshData GeometryGenerator::CreateCone() {
	MeshData meshData;
	int i = 0;

	float radius = 0.5f;
	int degrees = 2;
	float halfHeight = 0.5f;

	int baseTriangles = 360 / degrees;
	int triangles = baseTriangles * 2;
	int vertices = triangles * 3;
	meshData.Vertices.resize(vertices);
	meshData.Indices32.resize(vertices);

	for (int currentAngle = 0; currentAngle < 360; currentAngle += degrees) {

		float x1 = radius * cos(currentAngle * PI / 180);
		float x2 = radius * cos((currentAngle + degrees) * PI / 180);
		float z1 = radius * sin(currentAngle * PI / 180);
		float z2 = radius * sin((currentAngle + degrees) * PI / 180);

		/* DOWN Origin */
		meshData.Vertices[i++] = Vertex(0.0f, -halfHeight, 0.0f);
		/* A */
		meshData.Vertices[i++] = Vertex(x2,-halfHeight, z2);
		/* B */
		meshData.Vertices[i++] = Vertex(x1, -halfHeight, z1); 

		/* SIDE */
		/* A */
		meshData.Vertices[i++] = Vertex(0, halfHeight, 0);
		/* B */
		meshData.Vertices[i++] = Vertex(x1, -halfHeight, z1);
		/* C */
		meshData.Vertices[i++] = Vertex(x2, -halfHeight, z2); 

	}

	for (int j = 0; j < vertices; j++) {
		meshData.Indices32[j] = j;
	}

	return meshData;
}

GeometryGenerator::MeshData GeometryGenerator::CreateMySphere() {

	MeshData meshData;
	int i = 0;

	// SPHERE VERTICES VBO -----------------------------------------------------------
	float radius = 0.5f;
	int degrees = 10;
	float halfHeight = 0.5f;

	int baseTriangles = 360 / degrees;

	int triangles = baseTriangles * 2 /* TOP and BOTTOM */
		+ (baseTriangles - 2) * baseTriangles * 2; /* SIDES */
	int vertices = triangles * 3;
	meshData.Vertices.resize(vertices);
	meshData.Indices32.resize(vertices);

	for (int longAngle = 0; longAngle < 360; longAngle += degrees) {
		for (int latAngle = 0; latAngle < 360; latAngle += degrees) {
			float r1 = radius * sin(latAngle * PI / 180);
			float z1 = radius * cos(latAngle * PI / 180);
			float r1x1 = r1 * cos(longAngle * PI / 180);
			float r1x2 = r1 * cos((longAngle + degrees) * PI / 180);
			float r1z1 = r1 * sin(longAngle * PI / 180);
			float r1z2 = r1 * sin((longAngle + degrees) * PI / 180);

			float r2 = radius * sin((latAngle + degrees) * PI / 180);
			float z2 = radius * cos((latAngle + degrees) * PI / 180);
			float r2x1 = r2 * cos(longAngle * PI / 180);
			float r2x2 = r2 * cos((longAngle + degrees) * PI / 180);
			float r2z1 = r2 * sin(longAngle * PI / 180);
			float r2z2 = r2 * sin((longAngle + degrees) * PI / 180);

			if (latAngle == 0) {
				// TOP
				meshData.Vertices[i++] = Vertex(0.0f,radius, 0.0f);
				meshData.Vertices[i++] = Vertex(r2x1,z2,r2z1);
				meshData.Vertices[i++] = Vertex(r2x2,z2, r2z2);
			}
			else if (360 - latAngle <= degrees) {
				// BOTTOM
				meshData.Vertices[i++] = Vertex(0.0f,-radius,0.0f);
				meshData.Vertices[i++] = Vertex(r1x2, z1,r1z2);
				meshData.Vertices[i++] = Vertex(r1x1, z1, r1z1);
			}
			else {
				meshData.Vertices[i++] = Vertex(r1x2,z1,r1z2);
				meshData.Vertices[i++] = Vertex(r1x1,z1,r1z1);
				meshData.Vertices[i++] = Vertex(r2x1,z2,r2z1);
				meshData.Vertices[i++] = Vertex(r2x1,z2,r2z1);
				meshData.Vertices[i++] = Vertex(r2x2,z2,r2z2);
				meshData.Vertices[i++] = Vertex(r1x2,z1,r1z2);
			}
		}
	}

	for (int j = 0; j < vertices; j++) {
		meshData.Indices32[j] = j;
	}

	return meshData;
}
