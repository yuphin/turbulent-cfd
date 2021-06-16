#version 450
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
	float parity;
};

layout(binding = 5) readonly buffer CellType
{
	float cell_type[];
};

layout(binding = 6) buffer RS
{
	float rs[];
};

layout(binding = 7) readonly buffer P1
{
	float p1[];
};

layout(binding = 8) writeonly buffer P2
{
	float p2[];
};

float sor_helper(float p[5]) {
	float inv_dx2 = inv_dx * inv_dx;
	float inv_dy2 = inv_dy * inv_dy;
	float res = (p[0] + p[2]) * inv_dx2 + (p[3] + p[4]) * inv_dy2;
	return res;
}

void main() {
	uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    float is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 0) {
        return;
    }
	float p_stencil[5] = float[5](at(p1, i+1, j), at(p1, i, j), 
								  at(p1, i-1, j), at(p1, i, j+1), 
								  at(p1, i, j-1));
	float inv_dx2 = inv_dx * inv_dx;
	float inv_dy2 = inv_dy * inv_dy;
	float coeff = 1 / (2* (inv_dx2 + inv_dy2));
	at(p2, i, j) = at(p1, i, j) +
					   coeff * (sor_helper(p_stencil) - at(rs, i, j));
}