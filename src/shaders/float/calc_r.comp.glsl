#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "discretization.h"
#define at(AR, I, J) AR[imax * (J) + (I)]
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;


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

layout(binding = 14) buffer Residual
{
	float residual[];
};

void main() {
	uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
	float sum = 0.0;
    float is_fluid = at(cell_type, i, j);
	if((i < imax && j < jmax) && is_fluid == 1){
		float p_stencil[5] = float[5](at(p, i+1, j), at(p, i, j), 
								  at(p, i-1, j), at(p, i, j+1), 
								  at(p, i, j-1));
		float val = laplacian_5(p_stencil) - at(rs, i, j);
		sum = val * val;
		at(residual, i, j) = sum;
    }
	// sum = subgroupAdd(sum);
	// if (gl_SubgroupInvocationID == 0) {
    //     data[gl_SubgroupID] = sum;
	
    // }
	// barrier();
	
    // if (gl_SubgroupID == 0) {
    //     sum = data[gl_SubgroupInvocationID];
	// 	subgroupBarrier();
    //     sum = subgroupAdd(sum); 
    // }
	// if (gl_SubgroupInvocationID == 0) {
    //     r[gl_WorkGroupID.x] = sum;
    // }
}