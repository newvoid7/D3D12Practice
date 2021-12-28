//***************************************************************************************
// MainApp.cpp by Frank Luna (C) 2015 All Rights Reserved.
//***************************************************************************************

#include "../../Common/d3dApp.h"
#include "../../Common/MathHelper.h"
#include "../../Common/UploadBuffer.h"
#include "FrameResource.h"
#include "MeshData.h"

using Microsoft::WRL::ComPtr;
using namespace DirectX;
using namespace DirectX::PackedVector;

#pragma comment(lib, "d3dcompiler.lib")
#pragma comment(lib, "D3D12.lib")

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

enum class RenderLayer : int
{
    Opaque = 0,
    Count
};

class MainApp : public D3DApp
{
public:
    MainApp(HINSTANCE hInstance);
    MainApp(const MainApp& rhs) = delete;
    MainApp& operator=(const MainApp& rhs) = delete;
    ~MainApp();

    virtual bool Initialize()override;
    virtual LRESULT MsgProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)override;

private:
    virtual void OnResize()override;
    virtual void Update(const GameTimer& gt)override;
    virtual void Draw(const GameTimer& gt)override;

    virtual void OnMouseDown(WPARAM btnState, int x, int y)override;
    virtual void OnMouseUp(WPARAM btnState, int x, int y)override;
    virtual void OnMouseMove(WPARAM btnState, int x, int y)override;
    void OnMouseWheel(WPARAM wheel);
    void OnKeyDown(WPARAM key);
    void OnKeyUp(WPARAM key);
    
    void UpdateCamera(const GameTimer& gt);
    void UpdatePlane(const GameTimer& gt);
    void AnimateMaterials(const GameTimer& gt);
    void UpdateObjectCBs(const GameTimer& gt);
    void UpdateMaterialCBs(const GameTimer& gt);
    void UpdateMainPassCB(const GameTimer& gt);
    void UpdateWaves(const GameTimer& gt);

    void LoadTextures();
    void BuildRootSignature();
    void BuildDescriptorHeaps();
    void BuildShadersAndInputLayout();
    void BuildLandGeometry();
    void BuildWavesGeometry();
    void BuildPlaneGeometry();
    void BuildPSOs();
    void BuildFrameResources();
    void BuildMaterials();
    void BuildRenderItems();
    void DrawRenderItems(ID3D12GraphicsCommandList* cmdList, const std::vector<RenderItem*>& ritems);

    std::array<const CD3DX12_STATIC_SAMPLER_DESC, 6> GetStaticSamplers();

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

    std::vector<D3D12_INPUT_ELEMENT_DESC> mInputLayout;

    RenderItem* mWavesRitem = nullptr;

    // List of all the render items.
    std::vector<std::unique_ptr<RenderItem>> mAllRitems;

    // Render items divided by PSO.
    std::vector<RenderItem*> mRitemLayer[(int)RenderLayer::Count];

    std::unique_ptr<WaterMeshData> mWaves;

    PassConstants mMainPassCB;

private:
    // Interact related

    // 0: look down; 1: from the plane
    int mViewSet = 0;
    XMFLOAT3 mPlanePos = { 0.0f, 0.0f, 0.0f };
    XMFLOAT3 mPlaneAngle = { 0.0f, 0.0f, 0.0f }; // pitch, yaw, roll
    bool mPlaneForward = false, mPlaneBackward = false,
        mPlaneLeftward = false, mPlaneRightward = false,
        mPlaneLift = false, mPlaneSink = false;
    float mPlaneSpeed = 0.05f;
    XMFLOAT3 mEyePos = { 0.0f, 0.0f, 0.0f };
    //
    XMFLOAT3 mAtPos = { 0.0f, 0.0f, 0.0f };
    float mRadius = 100.0f;
    float mYaw = 0.0f;
    float mPitch = -XM_PIDIV4;
    //
    float mAzimuth = 0.0f;
    float mElevation = 0.0f;

    XMFLOAT4X4 mView = MathHelper::Identity4x4();
    XMFLOAT4X4 mProj = MathHelper::Identity4x4();

    //float mTheta = 1.5f * XM_PI;
    //float mPhi = XM_PIDIV2 - 0.1f;
    //float mRadius = 80.0f;

    POINT mLastMousePos;
};