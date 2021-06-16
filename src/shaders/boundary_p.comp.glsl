#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;

#define at(AR, I, J) AR[imax * (J) + (I)]

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

layout(binding = 21) buffer DT
{
	float dt;
};

layout(binding = 5) readonly buffer CellType
{
	float cell_type[];
};

layout(binding = 6) buffer RS
{
	float rs[];
};

layout(binding = 7) buffer P
{
	float p[];
};

layout(binding = 8) buffer R
{
	float r[];
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
	if(type == 0) {
		// Outlet
		at(p, i, j) = PI;
		//at(p, i, j) = imax * (j) + (i);
	} else  {
		// No-slip
		float p_stencil[8] = float[8](at(p, i+1, j), at(p, i - 1, j), 
								  at(p, i, j + 1), at(p, i, j-1), 
								  (at(p, i+1, j) + at(p, i, j + 1)) / 2,
								  (at(p, i+1, j) + at(p, i, j - 1)) / 2,
								  (at(p, i-1, j) + at(p, i, j + 1)) / 2,
								  (at(p, i-1, j) + at(p, i, j - 1)) / 2);
		at(p, i, j) = p_stencil[neighbors];
		//at(p, i, j) = imax * (j) + (i);//p_stencil[neighbors];
	}
}