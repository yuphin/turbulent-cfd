#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
layout(local_size_x = 1024, local_size_y = 1, local_size_z = 1) in;

layout(constant_id = 0) const uint MODE = 0;

layout(binding = 12) buffer V1 // d
{
	float v1[];
};

layout(binding = 13) buffer V2 // spmv
{
	float v2[];
};

layout(binding = 14) buffer Res // res
{
	float res[];
};

layout(binding = 15) buffer R // r
{
	float r[];
};

layout(binding = 16) buffer P // p
{
	float p[];
};

layout(binding = 17) buffer Z // z
{
	float z[];
};

void main() {
	uint idx = gl_GlobalInvocationID.x;
	// res = 0
	if(idx < size){
		float alpha =  res[0];
		switch(MODE){
			case 0:
				p[idx] = p[idx] + alpha * v1[idx];
				break;
			case 1:
				r[idx] = r[idx] - alpha * v2[idx];
				break;
			case 2:
				v1[idx] = r[idx] + alpha * v1[idx];
				break;
			case 3:
				v1[idx] = z[idx] + alpha * v1[idx];
				break;
		}
		
	} 
}