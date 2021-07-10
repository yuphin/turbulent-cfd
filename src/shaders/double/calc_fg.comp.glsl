#version 450
#extension GL_GOOGLE_include_directive: require
#include "discretization.h"
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;
layout(constant_id = 0) const uint CALC_TEMP = 0;
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

layout(binding = 5) readonly buffer CellType
{
	double cell_type[];
};

layout(binding = 21) buffer DT
{
	double dt;
};

layout(binding = 23) readonly buffer T
{
 	double t[];
};

void main() {
    uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    double is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 0){
        return;
    }
    // 5-point + 1 stencil for U and V
    double u_stencil[6] = double[6](at(u, i+1, j), at(u, i, j), at(u, i-1, j), at(u, i, j+1), at(u, i, j-1), at(u, i-1, j+1));
    double v_stencil[6] = double[6](at(v, i+1, j), at(v, i, j), at(v, i-1, j), at(v, i, j+1), at(v, i, j-1), at(v, i+1, j-1));

   // Calculate fluxes
   at(f, i, j) = at(u, i, j ) + dt * (nu * laplacian(u_stencil) - convection_u(u_stencil,v_stencil));
   at(g, i, j) = at(v, i, j ) + dt * (nu * laplacian(v_stencil) - convection_v(u_stencil,v_stencil));

   if(CALC_TEMP == 1){
       double term1 = at(t, i, j) + at(t, i+1, j);
       double term2 = at(t, i, j) + at(t, i, j + 1);
       at(f, i, j) -= beta * dt / 2 * (term1) * gx;
       at(g, i, j) -= beta * dt / 2 * (term2) * gy;
   }
}