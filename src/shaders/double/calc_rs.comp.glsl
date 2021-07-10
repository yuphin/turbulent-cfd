#version 450
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
#define at(AR, I, J) AR[imax * (J) + (I)]
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;

layout(std430, binding = 0) readonly buffer U
{
 	double u[];
};

layout(std430, binding = 1) readonly buffer V
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

double calculate_rs(double f[2], double g[2]) {
	double f_diff = inv_dx * ( f[0] - f[1]);
	double g_diff = inv_dy * ( g[0] - g[1]);
	return (f_diff + g_diff) * 1/ dt;
}


void main() {
    uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    double is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 0){
        return;
    }
	double f_diff[2] = double[2](at(f, i, j), at(f, i-1, j));
	double g_diff[2] = double[2](at(g, i, j), at(g, i, j - 1));
	at(rs, i, j) = calculate_rs(f_diff, g_diff);
	
}