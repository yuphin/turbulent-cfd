#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"

layout(local_size_x = 1024, local_size_y =1, local_size_z = 1) in;
layout(set = 1, binding = 21) buffer NutResidual
{
	double res[];
};

layout(binding = 31) buffer Counter
{
	int cnt;
};

shared double data[32];

void main(){
	uint idx = gl_GlobalInvocationID.x;
	double min_v = 0;
	int limit = cnt == 0 ? size : ((size + 1023) >> (10 * cnt));
	if(idx < limit){
		min_v = res[idx];
	}

	memoryBarrier();
	barrier();

	min_v = subgroupMin(min_v);
	if (gl_SubgroupInvocationID == 0) {
        data[gl_SubgroupID] = min_v;
    }
	barrier();

	if (gl_SubgroupID == 0) {
        min_v = data[gl_SubgroupInvocationID];
		subgroupBarrier();
        min_v = subgroupMin(min_v); 
    }
	if (gl_LocalInvocationID.x == 0) {
        res[gl_WorkGroupID.x] = min_v;
    }
	if(idx == 0){
		cnt = cnt + 1;
	}
}

