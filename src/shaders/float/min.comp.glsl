#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
layout(local_size_x = 1024, local_size_y =1, local_size_z = 1) in;
layout(set = 1, binding = 12) buffer NUT
{
	float NU_T[];
};

layout(binding = 14) buffer Residual //residual
{
	float res[];
};

shared float data[32];
void main(){
	uint idx = gl_GlobalInvocationID.x;
	float val = 0;
	if(idx < size){
		val = NU_T[idx];
	}
	val = subgroupMin(val);
	if (gl_SubgroupInvocationID == 0) {
        data[gl_SubgroupID] = val;
    }
	barrier();
	if (gl_SubgroupID == 0) {
    	val = data[gl_SubgroupInvocationID];
		subgroupBarrier();
        val = subgroupMax(val);
    }
	if (gl_LocalInvocationID.x == 0) {
        res[gl_WorkGroupID.x + 0] = val;
    }
}