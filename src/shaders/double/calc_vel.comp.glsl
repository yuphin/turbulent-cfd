#version 450
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
#define at(AR, I, J) AR[imax * (J) + (I)]
#define interpolate(A, i, j, i_offset, j_offset) (at(A, i, j) - at(A, i + i_offset, j + j_offset))/ 2
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;

layout(std430, binding = 0) buffer U
{
 	double u[];
};

layout(std430, binding = 1) buffer V
{
	double v[];
};

layout(std430, binding = 2) writeonly buffer F
{
	double f[];
};

layout(std430, binding = 3) writeonly buffer G
{
	double g[];
};

layout(binding = 21) buffer DT
{
	double dt;
};

layout(binding = 5) readonly buffer CellType
{
	double cell_type[];
};

layout(binding = 6) writeonly buffer RS
{
	double rs[];
};

layout(binding = 7) writeonly buffer P
{
	double p[];
};

double calculate_vel(double fg, double p[2], double inv_dxy) {
	
	return fg - dt * inv_dxy * (p[1] - p[0]);
}


void main() {
    uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    double is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 0){
        return;
    }
	double p_diff_u[2] = double[2](at(p, i, j), at(p, i+1, j));
	double p_diff_v[2] = double[2](at(p, i, j), at(p, i, j+1));
	at(u, i, j) = calculate_vel(at(f, i, j), p_diff_u, inv_dx);
	at(v, i, j) = calculate_vel(at(g, i, j), p_diff_v, inv_dy);
}