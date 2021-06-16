#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
#define at(AR, I, J) AR[imax * (J) + (I)]
#define interpolate(A, i, j, i_offset, j_offset) (at(A, i, j) - at(A, i + i_offset, j + j_offset))/ 2
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;

layout(std430, binding = 0) buffer U
{
 	float u[];
};

layout(std430, binding = 1) buffer V
{
	float v[];
};

layout(std430, binding = 2) buffer F
{
 	float f[];
};

layout(std430, binding = 3) buffer G
{
	float g[];
};

layout(binding = 5) readonly buffer CellType
{
	float cell_type[];
};

layout(binding = 9) buffer Neighborhood {
    uint neighborhood[];
};

void main() {
	uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    float is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 1) {
        return;
    }
	uint type = at(neighborhood, i, j) >> 8;
	uint neighbors = at(neighborhood, i, j) & 0xFF;
	if(neighbors == 0xFF) {
		return;
	}

	switch(neighbors){
			case 0: 
			{
				at(f, i, j) = at(u, i, j);
			}
			break;
			case 1: 
			{
				at(f, i-1, j) = at(u, i-1, j);
			}
			break;
			case 2: 
			{
				at(g, i, j) = at(v, i, j);
			}
			break;
			case 3: 
			{
				at(g, i, j-1) = at(v, i, j-1);
			}
			break;
	}
}