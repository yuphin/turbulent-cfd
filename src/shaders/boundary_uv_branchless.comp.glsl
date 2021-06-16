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

layout(binding = 22) buffer U
{
	float u[];
};

layout(binding = 23) buffer V
{
	float v[];
};

layout(binding = 24) buffer URHS
{
	float u_rhs[];
};

layout(binding = 25) buffer VRHS
{
	float v_rhs[];
};

layout(binding = 26) buffer URowStart
{
	int u_row_start[];
};

layout(binding = 27) buffer VRowStart
{
	int v_row_start[];
};

layout(binding = 28) buffer UColIdx
{
	int u_col_idx[];
};

layout(binding = 29) buffer VColIdx
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
		u_vec[row] = sum + u_rhs[row];

		sum = 0;
		for(int j = v_row_start[row]; j < v_row_start[row+1]; j++){
			sum += v[j] * v_vec[v_col_idx[j]];
		}
		v_vec[row] = sum + v_rhs[row];
	}
}
