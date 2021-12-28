#include <fstream>
#include <array>
#include <ppl.h>
#include <algorithm>
#include "../../Common/MathHelper.h"
#include "MeshData.h"

using std::wifstream;
using std::array;

wstring omit_starting_blank(const wstring& source)
{
	wstring::size_type i = 0;
	for (; i < source.size(); ++i) {
		if (source[i] != ' ' && source[i] != '\t' && source[i] != '\n') {
			break;
		}
	}
	return source.substr(i, source.size());
}

/// <summary>
/// Check if starts_with certain prefix, and set the post to following subwstring.
/// </summary>
/// <param name="source"></param>
/// <param name="start"></param>
/// <param name="post"></param>
/// <returns></returns>
bool starts_with(const wstring& source, const wstring& start, wstring& post)
{
	auto source_len = source.size(), start_len = start.size();
	if (source_len >= start_len) {
		post = source.substr(start_len, source.size());
		return source.substr(0, start_len) == start;
	}
	else {
		return false;
	}
}

bool ends_with(const wstring& source, const wstring& end)
{
	auto source_len = source.size(), end_len = end.size();
	if (source_len >= end_len) {
		return source.substr(source_len - end_len, source_len) == end;
	}
	else {
		return false;
	}
}

/// <summary>
/// https://www.zhihu.com/question/36642771/answer/865135551
/// </summary>
/// <param name="source"></param>
/// <param name="delimiters"></param>
/// <returns></returns>
vector<wstring> split(const wstring& source, const wstring& delimiters)
{
	vector<wstring> tokens;
	wstring::size_type last_pos = source.find_first_not_of(delimiters, 0);
	wstring::size_type pos = source.find_first_of(delimiters, last_pos);
	while (pos != wstring::npos || last_pos != wstring::npos) {
		tokens.push_back(source.substr(last_pos, pos - last_pos));
		last_pos = source.find_first_not_of(delimiters, pos);
		pos = source.find_first_of(delimiters, last_pos);
	}
	return tokens;
}

XMFLOAT2 parse_2_floats(const wstring& source, const wstring& delimiters)
{
	assert(!delimiters.empty());
	auto float_strs = split(source, delimiters);
	if (float_strs.size() < 2) {
		throw "Not enough floats.";
	};
	XMFLOAT2 ret(
		std::stof(float_strs[0]),
		std::stof(float_strs[1])
	);
	return ret;
}

XMFLOAT3 parse_3_floats(const wstring& source, const wstring& delimiters)
{
	assert(!delimiters.empty());
	auto float_strs = split(source, delimiters);
	if (float_strs.size() < 3) {
		throw "Not enough floats.";
	};
	XMFLOAT3 ret(
		std::stof(float_strs[0]),
		std::stof(float_strs[1]),
		std::stof(float_strs[2])
	);
	return ret;
}

array<array<std::uint16_t, 3>, 3> parse_3x3_uints(const wstring& source, const wstring& del0, const wstring& del1)
{
	assert(!del0.empty());
	auto uintx3_strs = split(source, del0);
	if (uintx3_strs.size() < 3) {
		throw "Not enough intx3.";
	}
	array<array<std::uint16_t, 3>, 3> ret;
	ret.fill(array<std::uint16_t, 3> { UINT16_MAX, UINT16_MAX, UINT16_MAX });
	for (int i = 0; i < 3; ++i) {
		assert(!del1.empty());
		auto s = split(uintx3_strs[i], del1);
		for (int j = 0; j < s.size(); ++j) {
			ret[i][j] = std::stoi(s[j]);
		}
	}
	return ret;
}

void MeshData::ReadMeshFile(const wstring& file_path)
{
	if (ends_with(file_path, L".stl")) {
		read_stl(file_path);
	}
	else if (ends_with(file_path, L".obj")) {
		read_obj(file_path);
	}
	else {

	}
}

void MeshData::compute_normals()
{
	vertices_normal = vector<XMFLOAT3>(VertexCount());
	XMVECTOR _0 = XMVectorSet(0.0f, 0.0f, 0.0f, 0.0f);
	vector<XMVECTOR> normal_help = vector<XMVECTOR>(VertexCount(), _0);
	for (size_t i = 0; i < indices.size() / 3; ++i) {
		auto p0 = XMLoadFloat3(&vertices_pos[indices[i * 3]]);
		auto p1 = XMLoadFloat3(&vertices_pos[indices[i * 3 + 1]]);
		auto p2 = XMLoadFloat3(&vertices_pos[indices[i * 3 + 2]]);
		auto normal = XMVector3Cross(p1 - p0, p2 - p0);
		normal_help[indices[i * 3]] += normal;
		normal_help[indices[i * 3 + 1]] += normal;
		normal_help[indices[i * 3 + 2]] += normal;
	}
	for (size_t i = 0; i < vertices_normal.size(); ++i) {
		XMStoreFloat3(&vertices_normal[i], XMVector3Normalize(normal_help[i]));
	}
}

/// <summary>
/// This will cover the pre-set vertices and indices
/// </summary>
/// <param name="file_path"></param>
void MeshData::read_stl(const wstring& file_path)
{
	wifstream ifs(file_path, std::ios::in);
	if (ifs.fail()) {
		return;
	}
	vector<XMFLOAT3> v_position;
	vector<std::uint16_t> v_index;
	// Record how many faces does a vertex belong, 
	// to compute the normal of this vertex
	XMFLOAT3 curr_normal;
	wstring line;
	while (!ifs.eof()) {
		std::getline(ifs, line);
		wstring post_str;

		if (starts_with(omit_starting_blank(line), L"facet normal ", post_str)) {
			curr_normal = parse_3_floats(post_str, L" ");
		}
		else if (starts_with(omit_starting_blank(line), L"vertex ", post_str)) {
			XMFLOAT3 p = parse_3_floats(post_str, L" ");
			auto find_iter = std::find_if(v_position.begin(), v_position.end(),
				[p](XMFLOAT3 it) {return it.x == p.x && it.y == p.y && it.z == p.z; });
			// If there doesn't exists the same vertex
			//  push back the new vertex
			if (find_iter == v_position.end()) {
				update_max_and_min(p);
				v_position.push_back(p);
				vertices_normal.push_back(curr_normal);
				v_index.push_back(v_position.size() - 1);
			}
			// If there already exits the same vertex
			//  add the count,
			//  add the sum of normal
			else {
				size_t index = find_iter - v_position.begin();
				v_index.push_back(index);
				vertices_normal[index].x += curr_normal.x;
				vertices_normal[index].y += curr_normal.y;
				vertices_normal[index].z += curr_normal.z;
			}
		}
	}
	vertices_pos = std::move(v_position);
	indices = std::move(v_index);
	// Read the normals
	for (size_t i = 0; i < vertices_normal.size(); ++i) {
		float norm = vertices_normal[i].x * vertices_normal[i].x
			+ vertices_normal[i].y * vertices_normal[i].y
			+ vertices_normal[i].z * vertices_normal[i].z;
		vertices_normal[i].x /= norm;
		vertices_normal[i].y /= norm;
		vertices_normal[i].z /= norm;
	}
	// Compute TexC
	vertices_texc.resize(VertexCount());
	for (size_t i = 0; i < VertexCount(); ++i) {
		vertices_texc[i].x = (vertices_pos[i].y + 1) / 2;
		vertices_texc[i].y = 1 - (vertices_pos[i].x + 1) / 2;
	}
}

void MeshData::read_obj(const wstring& file_path)
{
	wifstream ifs(file_path, std::ios::in);
	vector<std::uint16_t> n_index, t_index;
	vector<XMFLOAT3> v_normal;
	vector<XMFLOAT2> v_texc;
	wstring mtl_path;
	wstring line;
	while (!ifs.eof()) {
		std::getline(ifs, line);
		wstring post_str;
		if (starts_with(line, L"#", post_str)) {
			continue;
		}
		else {
			if (starts_with(line, L"v ", post_str)) {
				XMFLOAT3 p = parse_3_floats(post_str, L" ");
				update_max_and_min(p);
				vertices_pos.push_back(p);
			}
			else if (starts_with(line, L"vn ", post_str)) {
				XMFLOAT3 p = parse_3_floats(post_str, L" ");
				v_normal.push_back(p);
			}
			else if (starts_with(line, L"vt ", post_str)) {
				XMFLOAT2 p = parse_2_floats(post_str, L" ");
				v_texc.push_back(p);
			}
			else if (starts_with(line, L"f ", post_str)) {
				auto p = parse_3x3_uints(post_str, L" ", L"/");
				indices.push_back(p[0][0] - 1);
				indices.push_back(p[1][0] - 1);
				indices.push_back(p[2][0] - 1);
				t_index.push_back(p[0][1] - 1);
				t_index.push_back(p[1][1] - 1);
				t_index.push_back(p[2][1] - 1);
				n_index.push_back(p[0][2] - 1);
				n_index.push_back(p[1][2] - 1);
				n_index.push_back(p[2][2] - 1);
			}
			else if (starts_with(line, L"mtllib ", post_str)) {
				mtl_path = post_str;
			}
			else if (starts_with(line, L"usemtl ", post_str)) {

			}
			else if (starts_with(line, L"o ", post_str)
				|| starts_with(line, L"s ", post_str)) {

			}
			else {

			}
		}
	}
	vertices_normal.resize(VertexCount());
	vertices_texc.resize(VertexCount());
	concurrency::parallel_for(0, int(indices.size()), [this, v_normal, v_texc, n_index, t_index](int i) {
		vertices_normal[indices[i]] = n_index[i] != (UINT16_MAX - 1) ? 
			v_normal[n_index[i]] : XMFLOAT3(0.0f, 0.0f, 0.0f);
		vertices_texc[indices[i]] = t_index[i] != (UINT16_MAX-1) ? 
			v_texc[t_index[i]] : XMFLOAT2(0.0f, 0.0f);
	});
}

void MeshData::generate_grid(int num_x, int num_z)
{
	// count index
	auto index_count = [num_x, num_z](int i, int j) {return i * (num_z + 1) + j; };
	// vertex
	vertices_pos.resize((num_x + 1) * (num_z + 1));
	vertices_texc.resize((num_x + 1) * (num_z + 1));
	vertices_tangentu.resize((num_x + 1) * (num_z + 1));
	vertices_normal.resize((num_x + 1) * (num_z + 1));
	XMFLOAT3 min_position(-1, 0, -1), max_position(1, 0, 1);
	float dx = (max_position.x - min_position.x) / num_x;
	float dz = (max_position.z - min_position.z) / num_z;
	float x = min_position.x;
	int vi = 0;
	for (int i = 0; i < num_x + 1; ++i) {
		float z = min_position.z;
		for (int j = 0; j < num_z + 1; ++j) {
			vertices_pos[vi] = XMFLOAT3(x, 0, z);
			vertices_texc[vi] = XMFLOAT2(i * dx, j * dz);
			vertices_tangentu[vi] = XMFLOAT3(1.0f, 0.0f, 0.0f);
			vertices_tangentu[vi] = XMFLOAT3(0.0f, 1.0f, 0.0f);
			++vi;
			z += dz;
		}
		x += dx;
	}
	// index
	indices = vector<std::uint16_t>(num_x * num_z * 6);
	int ii = 0;
	for (int i = 0; i < num_x; ++i) {
		for (int j = 0; j < num_z; ++j) {
			indices[ii++] = index_count(i, j);
			indices[ii++] = index_count(i, j + 1);
			indices[ii++] = index_count(i + 1, j);
			indices[ii++] = index_count(i + 1, j);
			indices[ii++] = index_count(i, j + 1);
			indices[ii++] = index_count(i + 1, j + 1);
		}
	}
}

TerrainMeshData::TerrainMeshData(int num_x, int num_z)
	:MeshData(), cols(num_x), rows(num_z)
{
	generate_grid(num_x, num_z);
	// vertex magnitude
	for (auto& p : vertices_pos) {
		p.y = ComputeMagnitude(p.x, p.z);
	}
	// color
	vertices_color = vector<XMFLOAT4>(VertexCount());
	auto iter = vertices_pos.begin();
	for (auto& c : vertices_color) {
		auto y = (iter++)->y;
		if (y < -0.1f) {
			c = XMFLOAT4(1.0f, 0.96f, 0.62f, 1.0f);
		}
		else if (y < 0.1f) {
			c = XMFLOAT4(0.48f, 0.77f, 0.46f, 1.0f);
		}
		else if (y < 0.2f) {
			c = XMFLOAT4(0.1f, 0.48f, 0.19f, 1.0f);
		}
		else if (y < 0.25f) {
			c = XMFLOAT4(0.45f, 0.39f, 0.34f, 1.0f);
		}
		else {
			c = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);
		}
	}
	// normal
	compute_normals();
}

WaterMeshData::WaterMeshData(int num_x, int num_z, float time_step)
	:MeshData(), cols(num_x), rows(num_z), 
	time_step(time_step),
	spatial_step(2.0f / num_x)
{
	//float d = damping * time_step + 2.0f;
	//float e = (speed * speed) * (time_step * time_step) / (spatial_step * spatial_step);
	//k1 = (damping * time_step - 2.0f) / d;
	//k2 = (4.0f - 8.0f * e) / d;
	//k3 = (2.0f * e) / d;

	generate_grid(num_x, num_z);
	previous_pos = vertices_pos;
	vertices_pos = vertices_pos;
}

void WaterMeshData::UpdateWater(float dt)
{
	static float t = 0;

	// Accumulate time.
	t += dt;

	// Only update the simulation at the specified time step.
	if (t >= time_step)
	{
		// Only update interior points; we use zero boundary conditions.
		concurrency::parallel_for(1, rows, [this](int i)
			//for(int i = 1; i < mNumRows-1; ++i)
			{
				for (int j = 1; j < cols; ++j)
				{
					// After this update we will be discarding the old previous
					// buffer, so overwrite that buffer with the new update.
					// Note how we can do this inplace (read/write to same element) 
					// because we won't need prev_ij again and the assignment happens last.

					// Note j indexes x and i indexes z: h(x_j, z_i, t_k)
					// Moreover, our +z axis goes "down"; this is just to 
					// keep consistent with our row indices going down.

					previous_pos[IndexByRC(i, j)].y =
						k1 * previous_pos[IndexByRC(i, j)].y +
						k2 * vertices_pos[IndexByRC(i, j)].y +
						k3 * (vertices_pos[IndexByRC(i+1, j)].y +
							vertices_pos[IndexByRC(i-1, j)].y +
							vertices_pos[IndexByRC(i, j+1)].y +
							vertices_pos[IndexByRC(i, j-1)].y);
				}
			});

		// We just overwrote the previous buffer with the new data, so
		// this data needs to become the current solution and the old
		// current solution becomes the new previous solution.
		std::swap(previous_pos, vertices_pos);

		t = 0.0f; // reset time

		//
		// Compute normals using finite difference scheme.
		//
		concurrency::parallel_for(1, rows, [this](int i)
			//for(int i = 1; i < mNumRows - 1; ++i)
			{
				for (int j = 1; j < cols; ++j)
				{
					float l = vertices_pos[IndexByRC(i, j-1)].y;
					float r = vertices_pos[IndexByRC(i, j+1)].y;
					float t = vertices_pos[IndexByRC(i-1, j)].y;
					float b = vertices_pos[IndexByRC(i+1, j)].y;
					vertices_normal[IndexByRC(i, j)].x = -r + l;
					vertices_normal[IndexByRC(i, j)].y = 2.0f * spatial_step;
					vertices_normal[IndexByRC(i, j)].z = b - t;

					XMVECTOR n = XMVector3Normalize(XMLoadFloat3(&vertices_normal[IndexByRC(i, j)]));
					XMStoreFloat3(&vertices_normal[IndexByRC(i, j)], n);
				}
			});
	}
}

void WaterMeshData::Disturb(int i, int j, float magnitude)
{
	// Don't disturb boundaries.
	assert(i > 1 && i < rows - 1);
	assert(j > 1 && j < cols - 1);

	float halfMag = 0.5f * magnitude;

	// Disturb the ijth vertex height and its neighbors.
	vertices_pos[IndexByRC(i, j)].y += magnitude;
	vertices_pos[IndexByRC(i, j-1)].y += halfMag;
	vertices_pos[IndexByRC(i, j+1)].y += halfMag;
	vertices_pos[IndexByRC(i-1, j)].y += halfMag;
	vertices_pos[IndexByRC(i+1, j)].y += halfMag;
}

void WaterMeshData::RandomlyDisturb()
{
	int i = MathHelper::Rand(4, rows - 5);
	int j = MathHelper::Rand(4, cols - 5);
	float r = MathHelper::RandF(0.002f, 0.005f);
	Disturb(i, j, r);
}
