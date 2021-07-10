#version 450
#extension GL_GOOGLE_include_directive: require
#include "discretization.h"
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;
layout(constant_id = 0) const uint CALC_TEMP = 0;
layout(std430, binding = 0) readonly buffer U
{
 	float u[];
};

layout(std430, binding = 1) readonly buffer V
{
	float v[];
};

layout(std430, binding = 2) writeonly buffer F
{
	float f[];
};

layout(std430, binding = 3) writeonly buffer G
{
	float g[];
};

layout(binding = 5) readonly buffer CellType
{
	float cell_type[];
};

layout(binding = 21) buffer DT
{
	float dt;
};

layout(binding = 23) readonly buffer T
{
 	float t[];
};

void main() {
    uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    float is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 0){
        return;
    }
    // 5-point + 1 stencil for U and V
    float u_stencil[6] = float[6](at(u, i+1, j), at(u, i, j), at(u, i-1, j), at(u, i, j+1), at(u, i, j-1), at(u, i-1, j+1));
    float v_stencil[6] = float[6](at(v, i+1, j), at(v, i, j), at(v, i-1, j), at(v, i, j+1), at(v, i, j-1), at(v, i+1, j-1));

   // Calculate fluxes
   at(f, i, j) = at(u, i, j ) + dt * (nu * laplacian(u_stencil) - convection_u(u_stencil,v_stencil));
   at(g, i, j) = at(v, i, j ) + dt * (nu * laplacian(v_stencil) - convection_v(u_stencil,v_stencil));

   if(CALC_TEMP == 1){
       float term1 = at(t, i, j) + at(t, i+1, j);
       float term2 = at(t, i, j) + at(t, i, j + 1);
       at(f, i, j) -= beta * dt / 2 * (term1) * gx;
       at(g, i, j) -= beta * dt / 2 * (term2) * gy;
   }
}