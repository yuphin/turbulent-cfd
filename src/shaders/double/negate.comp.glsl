#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
layout(constant_id = 0) const uint STORE = 0;
layout(local_size_x = 1024, local_size_y = 1, local_size_z = 1) in;
layout(binding = 16) buffer P
{
	double p[];
};


void main(){
    uint idx = gl_GlobalInvocationID.x;
	if(idx < size){
        p[idx] = -p[idx];
    }
}

