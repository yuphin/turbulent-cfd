#version 450
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;

#define at(AR, I, J) AR[imax * (J) + (I)]


layout(std430, binding = 0) buffer U
{
 	double u[];
};

layout(std430, binding = 1) buffer V
{
	double v[];
};

layout(std430, binding = 2) buffer F
{
	double f[];
};

layout(std430, binding = 3) buffer G
{
	double g[];
};

layout(binding = 21) buffer DT
{
	double dt;
	double parity;
};

layout(binding = 5) readonly buffer CellType
{
	double cell_type[];
};

layout(binding = 6) buffer RS
{
	double rs[];
};

layout(binding = 7) readonly buffer P1
{
	double p1[];
};

layout(binding = 8) writeonly buffer P2
{
	double p2[];
};

double sor_helper(double p[5]) {
	double inv_dx2 = inv_dx * inv_dx;
	double inv_dy2 = inv_dy * inv_dy;
	double res = (p[0] + p[2]) * inv_dx2 + (p[3] + p[4]) * inv_dy2;
	return res;
}

void main() {
	uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    double is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 0) {
        return;
    }
	double p_stencil[5] = double[5](at(p1, i+1, j), at(p1, i, j), 
								  at(p1, i-1, j), at(p1, i, j+1), 
								  at(p1, i, j-1));
	double inv_dx2 = inv_dx * inv_dx;
	double inv_dy2 = inv_dy * inv_dy;
	double coeff = 1 / (2* (inv_dx2 + inv_dy2));
	at(p2, i, j) = at(p1, i, j) +
					   coeff * (sor_helper(p_stencil) - at(rs, i, j));
}