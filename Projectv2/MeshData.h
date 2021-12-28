#pragma once
#include <DirectXMath.h>
#include <DirectXColors.h>
#include <vector>
#include <string>
using namespace DirectX;
using std::vector;
using std::string;
using std::wstring;

class MeshData
{
protected:
    /// <summary>
    /// Not all vertex attributes will be used.
    /// If used, they should be the same size. (TODO: gaurantee that)
    /// If you want to pack some of them up, 
    /// do it just before writing to the buffer.
    /// </summary>
    vector<XMFLOAT3> vertices_pos;
    vector<XMFLOAT4> vertices_color;
    vector<XMFLOAT3> vertices_normal;
    vector<XMFLOAT3> vertices_tangentu;
    vector<XMFLOAT2> vertices_texc;
    vector<std::uint16_t> indices;
    vector<std::uint16_t> sub_indices;  // group of meshes
    vector<wstring> mtl_names;
public:
    MeshData() = default;
    void ReadMeshFile(const wstring& file_path);
public:
    int VertexCount() const { return vertices_pos.size(); }
    int IndexCount() const { return indices.size(); }
    void Border() {
        for (const auto& p : vertices_pos) {
            update_max_and_min(p);
        }
    };
    const vector<XMFLOAT3>& GetPositions() const { return vertices_pos; }
    const vector<XMFLOAT4>& GetColors() const { return vertices_color; }
    const vector<XMFLOAT3>& GetNormals() const { return vertices_normal; }
    const vector<XMFLOAT3>& GetTangentU() const { return vertices_tangentu; }
    const vector<XMFLOAT2>& GetTexC() const { return vertices_texc; }
    const vector<std::uint16_t>& GetIndices() const { return indices; }
protected:
    XMFLOAT3 Max_pos = { -FLT_MAX, -FLT_MAX, -FLT_MAX };
    XMFLOAT3 Min_pos = { FLT_MAX, FLT_MAX, FLT_MAX };
    void update_max_and_min(const XMFLOAT3& new_p) {
        Max_pos.x = XMMax(Max_pos.x, new_p.x);
        Max_pos.y = XMMax(Max_pos.y, new_p.y);
        Max_pos.z = XMMax(Max_pos.z, new_p.z);

        Min_pos.x = XMMin(Min_pos.x, new_p.x);
        Min_pos.y = XMMin(Min_pos.y, new_p.y);
        Min_pos.z = XMMin(Min_pos.z, new_p.z);
    }
    void compute_normals();
    void read_stl(const wstring& file_path);
    void read_obj(const wstring& file_path);
    void generate_grid(int num_x, int num_z);
};

class TerrainMeshData : public MeshData
{
    //TODO: load terrain from file
public:
    TerrainMeshData(int num_x = 200, int num_z = 200);
    float ComputeMagnitude(float x, float z) {
        float hill = 0.1f * (z * sinf(5.0f * x) + x * cosf(5.0f * z));
        float river = -0.3f;
        float distance_to_river = x - 0.2f * sinf(XM_PI * z);
        float w = expf(-10.0f * distance_to_river * distance_to_river);
        return w * river + (1 - w) * hill;
    }
private:
    int rows = 0;
    int cols = 0;
    int IndexByRC(int i, int j) { return i * (cols + 1) + j; }
};

class WaterMeshData : public MeshData
{
public:
    WaterMeshData(int num_x = 200, int num_z = 200,
        float time_step = 0.03f);
    void UpdateWater(float dt);
    void Disturb(int i, int j, float magnitude);
    void RandomlyDisturb();
private:
    int rows = 0;
    int cols = 0;

    // Simulation constants we can precompute.
    float k1 = -0.99402f;
    float k2 = 1.93659f;
    float k3 = 0.01435f;

    float time_step = 0.0f;
    float spatial_step = 0.0f;

    int IndexByRC(int i, int j) { return i * (cols + 1) + j; }

    vector<XMFLOAT3> previous_pos;
};