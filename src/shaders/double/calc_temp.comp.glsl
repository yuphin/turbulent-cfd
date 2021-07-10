#version 450
#extension GL_GOOGLE_include_directive: require
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;
#include "discretization.h"

layout(std430, binding = 0) readonly buffer U
{
 	double u[];
};

layout(std430, binding = 1) readonly buffer V
{
	double v[];
};

layout(binding = 5) readonly buffer CellType
{
	double cell_type[];
};

layout(binding = 21) buffer DT
{
	double dt;
};

layout(std430, binding = 22) readonly buffer T_old
{
 	double t_old[];
};

layout(std430, binding = 23) buffer T_new
{
 	double t_new[];
};

void main(){
    uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    double is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 0){
        return;
    }
    double u_stencil[2] = double[2](at(u, i-1, j), at(u, i, j));
    double v_stencil[2] = double[2](at(v, i, j - 1), at(v, i, j));

    double t_laplacian[5] = double[5](at(t_old, i+1, j), at(t_old, i, j), 
                    at(t_old, i-1, j), at(t_old, i, j+1), at(t_old, i, j-1));
                    
    at(t_new, i, j) = at(t_new, i, j) + dt * (alpha * laplacian_5(t_laplacian)
                    - convection_uT(u_stencil, t_laplacian)
                    - convection_vT(v_stencil, t_laplacian));
}
