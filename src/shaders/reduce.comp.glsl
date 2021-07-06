#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"

layout(local_size_x = 1024, local_size_y =1, local_size_z = 1) in;
layout(constant_id = 0) const uint STORE = 0;
layout(binding = 14) buffer Residual
{
	float res[];
};

layout(binding = 31) buffer Counter
{
	int cnt;
};

shared float data[32];

void main(){
	uint idx = gl_GlobalInvocationID.x;
	float sum = 0;
	int limit = cnt == 0 ? size : ((size + 1023) >> (10 * cnt));
	if(idx < limit){
		sum = res[idx];
	}

	memoryBarrier();
	barrier();

	sum = subgroupAdd(sum);
	if (gl_SubgroupInvocationID == 0) {
        data[gl_SubgroupID] = sum;
    }
	barrier();

	if (gl_SubgroupID == 0) {
        sum = data[gl_SubgroupInvocationID];
		subgroupBarrier();
        sum = subgroupAdd(sum); 
    }

	if (gl_LocalInvocationID.x == 0) {
        res[gl_WorkGroupID.x] = sum;
    }
	if(idx == 0){
		cnt = cnt + 1;
	}
}

