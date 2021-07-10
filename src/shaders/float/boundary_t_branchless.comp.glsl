#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
layout(local_size_x = 1024, local_size_y = 1, local_size_z = 1) in;

layout(binding = 23) buffer TVec
{
	float t_vec[];
};

layout(set = 1, binding = 8) buffer TBoundary
{
	float t[];
};

layout(set = 1, binding = 9) buffer TRHS
{
	float t_rhs[];
};

layout(set = 1, binding = 10) buffer TRowStart
{
	int t_row_start[];
};

layout(set = 1, binding = 11) buffer TColIdx
{
	int t_col_idx[];
};

void main() {
	uint row = gl_GlobalInvocationID.x;

	if(row < size) {
		float sum = 0;
		for(int j = t_row_start[row]; j < t_row_start[row+1]; j++){
			sum += t[j] * t_vec[t_col_idx[j]];
		}
		t_vec[row] = sum + 2 * t_rhs[row];
	}
}
