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
	float d[];
};


void main(){
	if(STORE == 0){
		res[0] =  d[0] / res[0];
	} else {
		d[1] = d[0]; // old
		d[0] = res[0]; // new
		res[0] = res[0] / d[1];
	}

}

