#include "MainApp.h"
#include "MeshData.h"
//#define DEBUG

const int gNumFrameResources = 3;

MainApp::MainApp(HINSTANCE hInstance)
	: D3DApp(hInstance)
{
	mMainWndCaption = L"Computer Graphics Project";
	mClientWidth = 1280;
	mClientHeight = 720;
}

MainApp::~MainApp()
{
	if (md3dDevice != nullptr)
		FlushCommandQueue();
}

bool MainApp::Initialize()
{
	if (!D3DApp::Initialize())
		return false;

	// Reset the command list to prep for initialization commands.
	ThrowIfFailed(mCommandList->Reset(mDirectCmdListAlloc.Get(), nullptr));

	// Get the increment size of a descriptor in this heap type.  This is hardware specific, 
	// so we have to query this information.
	mCbvSrvDescriptorSize = md3dDevice->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV);

	mWaves = std::make_unique<WaterMeshData>();

	LoadTextures();
	BuildRootSignature();
	BuildDescriptorHeaps();
	BuildShadersAndInputLayout();
	BuildLandGeometry();
	BuildWavesGeometry();
	BuildPlaneGeometry();
	BuildCAOGeometry();
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

LRESULT MainApp::MsgProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
	switch (msg)
	{
	case WM_KEYDOWN:
		OnKeyDown(wParam);
		return 0;
	case WM_KEYUP:
		OnKeyUp(wParam);
		return 0;
	case WM_MOUSEWHEEL:
		OnMouseWheel(wParam);
		return 0;
	default:
		return D3DApp::MsgProc(hwnd, msg, wParam, lParam);
	}
}

void MainApp::OnResize()
{
	D3DApp::OnResize();

	// The window resized, so update the aspect ratio and recompute the projection matrix.
	XMMATRIX P = XMMatrixPerspectiveFovLH(0.25f * MathHelper::Pi, AspectRatio(), 1.0f, 1000.0f);
	XMStoreFloat4x4(&mProj, P);
}

void MainApp::Update(const GameTimer& gt)
{
	UpdateCamera(gt);
	UpdatePlane(gt);

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

void MainApp::Draw(const GameTimer& gt)
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
	mCommandList->ClearRenderTargetView(CurrentBackBufferView(), Colors::LightSteelBlue, 0, nullptr);
	mCommandList->ClearDepthStencilView(DepthStencilView(), D3D12_CLEAR_FLAG_DEPTH | D3D12_CLEAR_FLAG_STENCIL, 1.0f, 0, 0, nullptr);

	// Specify the buffers we are going to render to.
	mCommandList->OMSetRenderTargets(1, &CurrentBackBufferView(), true, &DepthStencilView());

	ID3D12DescriptorHeap* descriptorHeaps[] = { mSrvDescriptorHeap.Get() };
	mCommandList->SetDescriptorHeaps(_countof(descriptorHeaps), descriptorHeaps);

	mCommandList->SetGraphicsRootSignature(mRootSignature.Get());

	auto passCB = mCurrFrameResource->PassCB->Resource();
	mCommandList->SetGraphicsRootConstantBufferView(2, passCB->GetGPUVirtualAddress());

	DrawRenderItems(mCommandList.Get(), mRitemLayer[(int)RenderLayer::Opaque]);

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

void MainApp::OnMouseDown(WPARAM btnState, int x, int y)
{
	mLastMousePos.x = x;
	mLastMousePos.y = y;

	SetCapture(mhMainWnd);
}

void MainApp::OnMouseUp(WPARAM btnState, int x, int y)
{
	ReleaseCapture();
}

void MainApp::OnMouseMove(WPARAM btnState, int x, int y)
{
	float dx = XMConvertToRadians(static_cast<float>(x - mLastMousePos.x));
	float dy = XMConvertToRadians(static_cast<float>(y - mLastMousePos.y));
	if (mViewMode == 0) {
		if ((btnState & MK_LBUTTON) != 0) {
			mAtPos.x += 1.5f * dy * cosf(mYaw) - 1.5f * dx * sinf(mYaw);
			mAtPos.z += 1.5f * dx * cosf(mYaw) + 1.5f * dy * sinf(mYaw);
		}
		else if ((btnState & MK_RBUTTON) != 0) {
		}
		else if ((btnState & MK_MBUTTON) != 0) {
			mYaw -= 0.2f * dx;
		}
	}
	else if (mViewMode == 1) {
		mAzimuth -= 360.0f * dy / mClientWidth;
		mElevation -= 90.0f * dx / mClientHeight;
		mElevation = MathHelper::Clamp(mElevation, -80.0f, 80.0f);
	}
	mLastMousePos.x = x;
	mLastMousePos.y = y;
}

void MainApp::OnMouseWheel(WPARAM wheel)
{
	short d = (short)HIWORD(wheel) / WHEEL_DELTA;
	if (mViewMode == 0) {
		mRadius -= 1.8f * d;
	}
	mRadius = MathHelper::Clamp(mRadius, 5.0f, 250.0f);
}

void MainApp::OnKeyDown(WPARAM key)
{
	switch (char(key)) {
	case 'W':
		mPlaneForward = true;
		break;
	case 'A':
		mPlaneLeftward = true;
		break;
	case 'S':
		mPlaneBackward = true;
		break;
	case 'D':
		mPlaneRightward = true;
		break;
	case 'Q':
		mPlaneLift = true;
		break;
	case 'E':
		mPlaneSink = true;
		break;
	case VK_SHIFT:
		mShiftPressed = true;
		break;
	case VK_TAB:
		mViewMode = 1 - mViewMode;
		break;
	case VK_SPACE:
		if (mViewMode == 0) {
			float dx = mEyePos.x - mAtPos.x, 
				dy = mEyePos.y - mAtPos.y, 
				dz = mEyePos.z - mAtPos.z;
			mAtPos.x = mPlanePos.x;
			mAtPos.y = mPlanePos.y - 10.0f;
			mAtPos.z = mPlanePos.z;
			mEyePos.x = mAtPos.x + dx;
			mEyePos.y = mAtPos.y + dy;
			mEyePos.z = mAtPos.z + dz;
		}
		break;
	case VK_UP:
		if (mViewMode == 0) {
			mAtPos.x += 1.0f * cosf(mYaw);
			mAtPos.z += 1.0f * sinf(mYaw);
		}
		break;
	case VK_DOWN:
		if (mViewMode == 0) {
			mAtPos.x -= 1.0f * cosf(mYaw);
			mAtPos.z -= 1.0f * sinf(mYaw);
		}
		break;
	case VK_LEFT:
		if (mViewMode == 0) {
			mAtPos.x -= 1.0f * sinf(mYaw);
			mAtPos.z += 1.0f * cosf(mYaw);
		}
		break;
	case VK_RIGHT:
		if (mViewMode == 0) {
			mAtPos.x += 1.0f * sinf(mYaw);
			mAtPos.z -= 1.0f * cosf(mYaw);
		}
		break;
	default:
		break;
	}
}

void MainApp::OnKeyUp(WPARAM key)
{
	switch (char(key))
	{
	case 'W':
		mPlaneForward = false;
		break;
	case 'A':
		mPlaneLeftward = false;
		break;
	case 'S':
		mPlaneBackward = false;
		break;
	case 'D':
		mPlaneRightward = false;
		break;
	case 'Q':
		mPlaneLift = false;
		break;
	case 'E':
		mPlaneSink = false;
		break;
	case VK_SHIFT:
		mShiftPressed = false;
		break;
	case VK_ESCAPE:
		PostQuitMessage(0);
		break;
	case VK_F2:
		//TODO
		// Set4xMsaaState(!m4xMsaaState);
		break;
	default:
		break;
	}
}

void MainApp::UpdateCamera(const GameTimer& gt)
{
	if (mViewMode == 0) {
		// Convert Spherical to Cartesian coordinates.
		mEyePos.x = mAtPos.x + mRadius * sinf(mPitch) * cosf(mYaw);
		mEyePos.z = mAtPos.z + mRadius * sinf(mPitch) * sinf(mYaw);
		mEyePos.y = mAtPos.y + mRadius * cosf(mPitch);
		//mEyePos = { 0.0f, mRadius, 0.0f };

		// Build the view matrix.
		XMVECTOR pos = XMVectorSet(mEyePos.x, mEyePos.y, mEyePos.z, 1.0f);
		XMVECTOR target = XMVectorSet(mAtPos.x, mAtPos.y, mAtPos.z, 1.0f);
		XMVECTOR up = XMVectorSet(0.0f, 1.0f, 0.0f, 0.0f);

		XMMATRIX view = XMMatrixLookAtLH(pos, target, up);
		XMStoreFloat4x4(&mView, view);
	}
	else if (mViewMode == 1) {
		mEyePos.x = mPlanePos.x;
		mEyePos.y = mPlanePos.y + 12.0f;
		mEyePos.z = mPlanePos.z;

		XMVECTOR pos = XMVectorSet(mEyePos.x, mEyePos.y, mEyePos.z, 1.0f);
		XMVECTOR target = XMVectorSet(
			cosf(mAzimuth) * cosf(mElevation) + mEyePos.x, 
			sinf(mAzimuth) * cosf(mElevation) + mEyePos.y, 
			sinf(mElevation) + mEyePos.z, 
			1.0f);
		XMVECTOR up = XMVectorSet(0.0f, 1.0f, 0.0f, 0.0f);
		
		XMMATRIX view = XMMatrixLookAtLH(pos, target, up);
		XMStoreFloat4x4(&mView, view);
	}
}

void MainApp::UpdatePlane(const GameTimer& gt)
{
	// Use a timer to control the update
	if (mPlaneForward || mPlaneBackward
		|| mPlaneLeftward || mPlaneRightward
		|| mPlaneLift || mPlaneSink) {
		if (!mShiftPressed) {
			if (mPlaneForward) {
				mPlanePos.x += mPlaneSpeed;
			}
			if (mPlaneBackward) {
				mPlanePos.x -= mPlaneSpeed;
			}
			if (mPlaneLeftward) {
				mPlanePos.z += mPlaneSpeed;
			}
			if (mPlaneRightward) {
				mPlanePos.z -= mPlaneSpeed;
			}
			if (mPlaneLift) {
				mPlanePos.y += mPlaneSpeed;
			}
			if (mPlaneSink) {
				mPlanePos.y -= mPlaneSpeed;
			}
		}
		else {
			if (mPlaneForward) {
				mPlaneAngle.y += mPlaneRotateSpeed;
			}
			if (mPlaneBackward) {
				mPlaneAngle.y -= mPlaneRotateSpeed;
			}
			if (mPlaneLeftward) {
				mPlaneAngle.z -= mPlaneRotateSpeed;
			}
			if (mPlaneRightward) {
				mPlaneAngle.z += mPlaneRotateSpeed;
			}
			if (mPlaneLift) {
				mPlaneAngle.x += mPlaneRotateSpeed;
			}
			if (mPlaneSink) {
				mPlaneAngle.x -= mPlaneRotateSpeed;
			}
		}
		// Important!
		mAllRitems[2]->NumFramesDirty = gNumFrameResources;
	}
}

void MainApp::AnimateMaterials(const GameTimer& gt)
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

void MainApp::UpdateObjectCBs(const GameTimer& gt)
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

			if (e->Geo->Name == "planeGeo") {
				world = XMMatrixMultiply(
					XMMatrixRotationRollPitchYaw(mPlaneAngle.x, mPlaneAngle.y, mPlaneAngle.z),
					XMMatrixMultiply(
						world,
						XMMatrixTranslation(mPlanePos.x, mPlanePos.y, mPlanePos.z)
					)
				);
			}

			ObjectConstants objConstants;
			XMStoreFloat4x4(&objConstants.World, XMMatrixTranspose(world));
			XMStoreFloat4x4(&objConstants.TexTransform, XMMatrixTranspose(texTransform));

			currObjectCB->CopyData(e->ObjCBIndex, objConstants);

			// Next FrameResource need to be updated too.
			e->NumFramesDirty--;
		}
	}
}

void MainApp::UpdateMaterialCBs(const GameTimer& gt)
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

void MainApp::UpdateMainPassCB(const GameTimer& gt)
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
	mMainPassCB.AmbientLight = { 0.25f, 0.25f, 0.35f, 1.0f };

	XMFLOAT3 sunDirection = {
		0.707106f * sinf(0.8f * gt.TotalTime()),
		-0.707106f * cosf(0.8f * gt.TotalTime()),
		0.707106f
	};

	//mMainPassCB.Lights[0].Direction = { 0.57735f, -0.57735f, 0.57735f };
	mMainPassCB.Lights[0].Direction = { sunDirection.x, sunDirection.y, sunDirection.z };
	mMainPassCB.Lights[0].Strength = { 0.9f, 0.9f, 0.9f };
	//mMainPassCB.Lights[1].Direction = { -0.57735f, -0.57735f, 0.57735f };
	mMainPassCB.Lights[1].Direction = { -sunDirection.x, -sunDirection.y, -sunDirection.z };
	mMainPassCB.Lights[1].Strength = { 0.2f, 0.2f, 0.8f };
	mMainPassCB.Lights[2].Direction = { 0.0f, -0.707f, -0.707f };
	mMainPassCB.Lights[2].Strength = { 0.2f, 0.2f, 0.2f };

	auto currPassCB = mCurrFrameResource->PassCB.get();
	currPassCB->CopyData(0, mMainPassCB);
}

void MainApp::UpdateWaves(const GameTimer& gt)
{
	// Every quarter second, generate a random wave.
	static float t_base = 0.0f;
	if ((mTimer.TotalTime() - t_base) >= 0.25f)
	{
		t_base += 0.25f;
		mWaves->RandomlyDisturb();
	}

	// Update the wave simulation.
	mWaves->UpdateWater(gt.DeltaTime());

	// Update the wave vertex buffer with the new solution.
	auto currWavesVB = mCurrFrameResource->WavesVB.get();
	for (int i = 0; i < mWaves->VertexCount(); ++i)
	{
		Vertex v;

		v.Pos = mWaves->GetPositions()[i];
		v.Normal = mWaves->GetNormals()[i];
		v.TexC = mWaves->GetTexC()[i];

		currWavesVB->CopyData(i, v);
	}

	// Set the dynamic VB of the wave renderitem to the current frame VB.
	mWavesRitem->Geo->VertexBufferGPU = currWavesVB->Resource();
}

void MainApp::LoadTextures()
{
	auto grassTex = std::make_unique<Texture>();
	grassTex->Name = "grassTex";
#ifdef DEBUG
	grassTex->Filename = L"../../Textures/grass.dds";
#endif // DEBUG
#ifndef DEBUG
	grassTex->Filename = L"Data/grass.dds";
#endif // !DEBUG
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(md3dDevice.Get(),
		mCommandList.Get(), grassTex->Filename.c_str(),
		grassTex->Resource, grassTex->UploadHeap));

	auto waterTex = std::make_unique<Texture>();
	waterTex->Name = "waterTex";
#ifdef DEBUG
	waterTex->Filename = L"../../Textures/water1.dds";
#endif // DEBUG
#ifndef DEBUG
	waterTex->Filename = L"Data/water1.dds";
#endif // !DEBUG
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(md3dDevice.Get(),
		mCommandList.Get(), waterTex->Filename.c_str(),
		waterTex->Resource, waterTex->UploadHeap));

	auto planeTex = std::make_unique<Texture>();
	planeTex->Name = "planeTex";
#ifdef DEBUG
	planeTex->Filename = L"../../Textures/plane.dds";
#endif // DEBUG
#ifndef DEBUG
	planeTex->Filename = L"Data/plane.dds";
#endif // !DEBUG
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(md3dDevice.Get(),
		mCommandList.Get(), planeTex->Filename.c_str(),
		planeTex->Resource, planeTex->UploadHeap));

	auto caoTex = std::make_unique<Texture>();
	caoTex->Name = "caoTex";
#ifdef DEBUG
	caoTex->Filename = L"../../Textures/other.dds";
#endif // DEBUG
#ifndef DEBUG
	caoTex->Filename = L"Data/other.dds";
#endif // !DEBUG
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(md3dDevice.Get(),
		mCommandList.Get(), caoTex->Filename.c_str(),
		caoTex->Resource, caoTex->UploadHeap));

	mTextures[grassTex->Name] = std::move(grassTex);
	mTextures[waterTex->Name] = std::move(waterTex);
	mTextures[planeTex->Name] = std::move(planeTex);
	mTextures[caoTex->Name] = std::move(caoTex);
}

void MainApp::BuildRootSignature()
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

void MainApp::BuildDescriptorHeaps()
{
	//
	// Create the SRV heap.
	//
	D3D12_DESCRIPTOR_HEAP_DESC srvHeapDesc = {};
	srvHeapDesc.NumDescriptors = 4;
	srvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
	srvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;
	ThrowIfFailed(md3dDevice->CreateDescriptorHeap(&srvHeapDesc, IID_PPV_ARGS(&mSrvDescriptorHeap)));

	//
	// Fill out the heap with actual descriptors.
	//
	CD3DX12_CPU_DESCRIPTOR_HANDLE hDescriptor(mSrvDescriptorHeap->GetCPUDescriptorHandleForHeapStart());

	auto grassTex = mTextures["grassTex"]->Resource;
	auto waterTex = mTextures["waterTex"]->Resource;
	auto planeTex = mTextures["planeTex"]->Resource;
	auto caoTex = mTextures["caoTex"]->Resource;

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

	srvDesc.Format = planeTex->GetDesc().Format;
	md3dDevice->CreateShaderResourceView(planeTex.Get(), &srvDesc, hDescriptor);

	// next descriptor
	hDescriptor.Offset(1, mCbvSrvDescriptorSize);

	srvDesc.Format = caoTex->GetDesc().Format;
	md3dDevice->CreateShaderResourceView(caoTex.Get(), &srvDesc, hDescriptor);

}

void MainApp::BuildShadersAndInputLayout()
{
#ifdef DEBUG
	mShaders["standardVS"] = d3dUtil::CompileShader(L"Shaders\\Default.hlsl", nullptr, "VS", "vs_5_0");
	mShaders["opaquePS"] = d3dUtil::CompileShader(L"Shaders\\Default.hlsl", nullptr, "PS", "ps_5_0");
#endif // DEBUG
#ifndef DEBUG
	mShaders["standardVS"] = d3dUtil::CompileShader(L"Data\\Default.hlsl", nullptr, "VS", "vs_5_0");
	mShaders["opaquePS"] = d3dUtil::CompileShader(L"Data\\Default.hlsl", nullptr, "PS", "ps_5_0");
#endif // !DEBUG

	mInputLayout =
	{
		{ "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
		{ "NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
		{ "TEXCOORD", 0, DXGI_FORMAT_R32G32_FLOAT, 0, 24, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
	};
}

void MainApp::BuildLandGeometry()
{
	TerrainMeshData land;
	vector<Vertex> vertices(land.VertexCount());
	for (size_t i = 0; i < land.VertexCount(); ++i) {
		vertices[i].Pos = land.GetPositions()[i];
		vertices[i].Normal = land.GetNormals()[i];
		vertices[i].TexC = land.GetTexC()[i];
	}
	std::vector<std::uint16_t> indices = land.GetIndices();

	const UINT vbByteSize = (UINT)vertices.size() * sizeof(Vertex);
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

void MainApp::BuildWavesGeometry()
{
	UINT vbByteSize = mWaves->VertexCount() * sizeof(Vertex);
	//UINT ibByteSize = (UINT)indices.size() * sizeof(std::uint16_t);
	UINT ibByteSize = (UINT)mWaves->IndexCount() * sizeof(std::uint16_t);

	auto geo = std::make_unique<MeshGeometry>();
	geo->Name = "waterGeo";

	// Set dynamically.
	geo->VertexBufferCPU = nullptr;
	geo->VertexBufferGPU = nullptr;

	ThrowIfFailed(D3DCreateBlob(ibByteSize, &geo->IndexBufferCPU));
	//CopyMemory(geo->IndexBufferCPU->GetBufferPointer(), indices.data(), ibByteSize);
	CopyMemory(geo->IndexBufferCPU->GetBufferPointer(), mWaves->GetIndices().data(), ibByteSize);

	geo->IndexBufferGPU = d3dUtil::CreateDefaultBuffer(md3dDevice.Get(),
		//mCommandList.Get(), indices.data(), ibByteSize, geo->IndexBufferUploader);
		mCommandList.Get(), mWaves->GetIndices().data(), ibByteSize, geo->IndexBufferUploader);

	geo->VertexByteStride = sizeof(Vertex);
	geo->VertexBufferByteSize = vbByteSize;
	geo->IndexFormat = DXGI_FORMAT_R16_UINT;
	geo->IndexBufferByteSize = ibByteSize;

	SubmeshGeometry submesh;
	//submesh.IndexCount = (UINT)indices.size();
	submesh.IndexCount = (UINT)mWaves->IndexCount();
	submesh.StartIndexLocation = 0;
	submesh.BaseVertexLocation = 0;

	geo->DrawArgs["grid"] = submesh;

	mGeometries["waterGeo"] = std::move(geo);
}

void MainApp::BuildPlaneGeometry()
{
	MeshData plane;
#ifdef DEBUG
	plane.ReadMeshFile(L"..\\..\\..\\model\\f117.stl");
#endif // DEBUG
#ifndef DEBUG
	plane.ReadMeshFile(L"Data\\f117.stl");
#endif // !DEBUG
	std::vector<Vertex> vertices(plane.VertexCount());
	std::vector<std::uint16_t> indices = plane.GetIndices();
	size_t i = 0;
	for (auto& v : vertices) {
		v.Pos = plane.GetPositions()[i];
		v.Normal = plane.GetNormals()[i];
		v.TexC = plane.GetTexC()[i];
		++i;
	}

	const UINT vbByteSize = (UINT)vertices.size() * sizeof(Vertex);
	const UINT ibByteSize = (UINT)indices.size() * sizeof(std::uint16_t);

	auto geo = std::make_unique<MeshGeometry>();
	geo->Name = "planeGeo";

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

	geo->DrawArgs["plane"] = submesh;

	mGeometries["planeGeo"] = std::move(geo);
}

void MainApp::BuildCAOGeometry()
{
	MeshData cao;
#ifdef DEBUG
	cao.ReadMeshFile(L"..\\..\\..\\model\\blender_grass.obj");
#endif // DEBUG
#ifndef DEBUG
	cao.ReadMeshFile(L"Data\\blender_grass.obj");
#endif // !DEBUG
	std::vector<Vertex> vertices(cao.VertexCount());
	std::vector<std::uint16_t> indices = cao.GetIndices();
	size_t i = 0;
	for (auto& v : vertices) {
		v.Pos = cao.GetPositions()[i];
		v.Normal = cao.GetNormals()[i];
		v.TexC = cao.GetTexC()[i];
		++i;
	}

	const UINT vbByteSize = (UINT)vertices.size() * sizeof(Vertex);
	const UINT ibByteSize = (UINT)indices.size() * sizeof(std::uint16_t);

	auto geo = std::make_unique<MeshGeometry>();
	geo->Name = "caoGeo";

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

	geo->DrawArgs["cao"] = submesh;

	mGeometries["caoGeo"] = std::move(geo);
}

void MainApp::BuildPSOs()
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
	ThrowIfFailed(md3dDevice->CreateGraphicsPipelineState(&opaquePsoDesc, IID_PPV_ARGS(&mPSOs["opaque"])));
}

void MainApp::BuildFrameResources()
{
	for (int i = 0; i < gNumFrameResources; ++i)
	{
		mFrameResources.push_back(std::make_unique<FrameResource>(md3dDevice.Get(),
			1, (UINT)mAllRitems.size(), (UINT)mMaterials.size(), mWaves->VertexCount()));
	}
}

void MainApp::BuildMaterials()
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
	water->DiffuseAlbedo = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);
	water->FresnelR0 = XMFLOAT3(0.2f, 0.2f, 0.2f);
	water->Roughness = 0.0f;

	auto plane = std::make_unique<Material>();
	plane->Name = "plane";
	plane->MatCBIndex = 2;
	plane->DiffuseSrvHeapIndex = 2;
	plane->DiffuseAlbedo = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);
	plane->FresnelR0 = XMFLOAT3(0.1f, 0.1f, 0.1f);
	plane->Roughness = 0.25f;

	auto cao = std::make_unique<Material>();
	cao->Name = "cao";
	cao->MatCBIndex = 3;
	cao->DiffuseSrvHeapIndex = 3;
	cao->DiffuseAlbedo = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);
	cao->FresnelR0 = XMFLOAT3(0.1f, 0.1f, 0.1f);
	cao->Roughness = 0.05f;

	mMaterials["grass"] = std::move(grass);
	mMaterials["water"] = std::move(water);
	mMaterials["plane"] = std::move(plane);
	mMaterials["cao"] = std::move(cao);
}

void MainApp::BuildRenderItems()
{
	auto wavesRitem = std::make_unique<RenderItem>();
	XMStoreFloat4x4(&wavesRitem->World,
		XMMatrixMultiply(
			XMMatrixScaling(60.0f, 60.0f, 60.0f),
			XMMatrixTranslation(0.0f, 0.0f, 0.0f)
		)
	);
	XMStoreFloat4x4(&wavesRitem->TexTransform, XMMatrixScaling(5.0f, 5.0f, 1.0f));
	wavesRitem->ObjCBIndex = 0;
	wavesRitem->Mat = mMaterials["water"].get();
	wavesRitem->Geo = mGeometries["waterGeo"].get();
	wavesRitem->PrimitiveType = D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST;
	wavesRitem->IndexCount = wavesRitem->Geo->DrawArgs["grid"].IndexCount;
	wavesRitem->StartIndexLocation = wavesRitem->Geo->DrawArgs["grid"].StartIndexLocation;
	wavesRitem->BaseVertexLocation = wavesRitem->Geo->DrawArgs["grid"].BaseVertexLocation;
	mWavesRitem = wavesRitem.get();
	mRitemLayer[(int)RenderLayer::Opaque].push_back(wavesRitem.get());
	mAllRitems.push_back(std::move(wavesRitem));

	auto terrainRitem = std::make_unique<RenderItem>();
	XMStoreFloat4x4(&terrainRitem->World,
		XMMatrixMultiply(
			XMMatrixScaling(60.0f, 100.0f, 60.0f),
			XMMatrixTranslation(0.0, 20.0f, 0.0f)
		)
	);
	XMStoreFloat4x4(&terrainRitem->TexTransform, XMMatrixScaling(5.0f, 5.0f, 1.0f));
	terrainRitem->ObjCBIndex = 1;
	terrainRitem->Mat = mMaterials["grass"].get();
	terrainRitem->Geo = mGeometries["landGeo"].get();
	terrainRitem->PrimitiveType = D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST;
	terrainRitem->IndexCount = terrainRitem->Geo->DrawArgs["grid"].IndexCount;
	terrainRitem->StartIndexLocation = terrainRitem->Geo->DrawArgs["grid"].StartIndexLocation;
	terrainRitem->BaseVertexLocation = terrainRitem->Geo->DrawArgs["grid"].BaseVertexLocation;
	mRitemLayer[(int)RenderLayer::Opaque].push_back(terrainRitem.get());
	mAllRitems.push_back(std::move(terrainRitem));

	auto planeRitem = std::make_unique<RenderItem>();
	XMStoreFloat4x4(&planeRitem->World,
		XMMatrixMultiply(
			XMMatrixMultiply(XMMatrixRotationX(-XM_PIDIV2), XMMatrixScaling(10.0f, 10.0f, 10.0f)),
			XMMatrixTranslation(0.0f, 10.0f, 0.0f)
		)
	);
	planeRitem->ObjCBIndex = 2;
	planeRitem->Mat = mMaterials["plane"].get();
	planeRitem->Geo = mGeometries["planeGeo"].get();
	planeRitem->PrimitiveType = D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST;
	planeRitem->IndexCount = planeRitem->Geo->DrawArgs["plane"].IndexCount;
	planeRitem->StartIndexLocation = planeRitem->Geo->DrawArgs["plane"].StartIndexLocation;
	planeRitem->BaseVertexLocation = planeRitem->Geo->DrawArgs["plane"].BaseVertexLocation;
	mRitemLayer[(int)RenderLayer::Opaque].push_back(planeRitem.get());
	mAllRitems.push_back(std::move(planeRitem));

	auto caoRitem = std::make_unique<RenderItem>();
	XMStoreFloat4x4(&caoRitem->World,
		XMMatrixMultiply(
			XMMatrixScaling(1.0f, 1.0f, 1.0f),
			XMMatrixTranslation(30.0f, 10.0f, -45.0f)
		)
	);
	caoRitem->ObjCBIndex = 3;
	caoRitem->Mat = mMaterials["cao"].get();
	caoRitem->Geo = mGeometries["caoGeo"].get();
	caoRitem->PrimitiveType = D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST;
	caoRitem->IndexCount = caoRitem->Geo->DrawArgs["cao"].IndexCount;
	caoRitem->StartIndexLocation = caoRitem->Geo->DrawArgs["cao"].StartIndexLocation;
	caoRitem->BaseVertexLocation = caoRitem->Geo->DrawArgs["cao"].BaseVertexLocation;
	mRitemLayer[(int)RenderLayer::Opaque].push_back(caoRitem.get());
	mAllRitems.push_back(std::move(caoRitem));
}

void MainApp::DrawRenderItems(ID3D12GraphicsCommandList* cmdList, const std::vector<RenderItem*>& ritems)
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

		CD3DX12_GPU_DESCRIPTOR_HANDLE tex(mSrvDescriptorHeap->GetGPUDescriptorHandleForHeapStart());
		tex.Offset(ri->Mat->DiffuseSrvHeapIndex, mCbvSrvDescriptorSize);

		D3D12_GPU_VIRTUAL_ADDRESS objCBAddress = objectCB->GetGPUVirtualAddress() + ri->ObjCBIndex * objCBByteSize;
		D3D12_GPU_VIRTUAL_ADDRESS matCBAddress = matCB->GetGPUVirtualAddress() + ri->Mat->MatCBIndex * matCBByteSize;

		cmdList->SetGraphicsRootDescriptorTable(0, tex);
		cmdList->SetGraphicsRootConstantBufferView(1, objCBAddress);
		cmdList->SetGraphicsRootConstantBufferView(3, matCBAddress);

		cmdList->DrawIndexedInstanced(ri->IndexCount, 1, ri->StartIndexLocation, ri->BaseVertexLocation, 0);
	}
}

std::array<const CD3DX12_STATIC_SAMPLER_DESC, 6> MainApp::GetStaticSamplers()
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
