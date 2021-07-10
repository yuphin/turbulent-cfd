#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
layout(constant_id = 0) const uint STORE = 0;
layout(local_size_x = 1, local_size_y =1, local_size_z = 1) in;

layout(binding = 14) buffer Residual
{
	float res[];
};

layout(binding = 30) buffer Deltas
{
	float deltas[];
};


void main(){
	if(STORE == 0){
		res[0] =  deltas[0] / res[0];
	} else if(STORE == 1){
		deltas[1] = deltas[0]; // old
		deltas[0] = res[0]; // new
		res[0] = res[0] / deltas[1];
	} else if(STORE == 2){
		res[0] = sqrt(res[0] / num_fluid_cells);
	}

}

