#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
layout(local_size_x = 1024, local_size_y = 1, local_size_z = 1) in;


layout(binding = 0) buffer U
{
 	float u[];
};

layout(binding = 1) buffer V
{
	float v[];
};

layout(binding = 2) buffer F
{
 	float f[];
};

layout(binding = 3) buffer G
{
	float g[];
};

layout(binding = 5) readonly buffer CellType
{
	float cell_type[];
};


layout(set = 1, binding = 12) buffer FIdx
{
	int f_idx[];
};

layout(set = 1, binding = 13) buffer GIdx
{
	int g_idx[];
};
void main() {
	uint id = gl_GlobalInvocationID.x;
	float is_fluid = cell_type[id];
	if(id < size && is_fluid != 1) {
		int uidx = f_idx[id];
		int vidx = g_idx[id];
		if(uidx != -1){
			f[uidx] = u[uidx];
		}
		if(vidx != -1){
			g[vidx] = v[vidx];
		}
		
	}
}
