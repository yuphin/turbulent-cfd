#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
layout(local_size_x = 1024, local_size_y = 1, local_size_z = 1) in;

layout(binding = 0) buffer UVec
{
	float u_vec[];
};

layout(binding = 1) buffer Vvec
{
	float v_vec[];
};

layout(set = 1, binding = 0) buffer U
{
	float u[];
};

layout(set = 1, binding = 1) buffer V
{
	float v[];
};

layout(set = 1, binding = 2) buffer URHS
{
	float u_rhs[];
};

layout(set = 1, binding = 3) buffer VRHS
{
	float v_rhs[];
};

layout(set = 1, binding = 4) buffer URowStart
{
	int u_row_start[];
};

layout(set = 1, binding = 5) buffer VRowStart
{
	int v_row_start[];
};

layout(set = 1, binding = 6) buffer UColIdx
{
	int u_col_idx[];
};

layout(set = 1, binding = 7) buffer VColIdx
{
	int v_col_idx[];
};

void main() {
	uint row = gl_GlobalInvocationID.x;

	if(row < size) {
		float sum = 0;
		for(int j = u_row_start[row]; j < u_row_start[row+1]; j++){
			sum += u[j] * u_vec[u_col_idx[j]];
		}
		u_vec[row] = sum + 2 * u_rhs[row];

		sum = 0;
		for(int j = v_row_start[row]; j < v_row_start[row+1]; j++){
			sum += v[j] * v_vec[v_col_idx[j]];
		}
		v_vec[row] = sum + 2 * v_rhs[row];
	}
}
