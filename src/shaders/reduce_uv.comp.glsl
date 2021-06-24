#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"

layout(local_size_x = 1024, local_size_y =1, local_size_z = 1) in;
layout(constant_id = 0) const uint UV = 0;
layout(binding = 14) buffer Residual //residual
{
	float res[];
};

layout(binding = 31) buffer Counter
{
	int cnt;
};

shared float data_max[32];
void main(){
	uint idx = gl_GlobalInvocationID.x;
	float max_value = 0;
	int limit = ((size + 1023) >> (10 * cnt)); 
	if(idx < limit){
		if(UV == 0){
		    max_value = res[2*idx + 0];
		} else if(UV == 1){
            max_value = res[2*idx + 1];
		}
	}
	memoryBarrier();
	barrier();
	
	max_value = subgroupMax(max_value);

	if (gl_SubgroupInvocationID == 0) {
        data_max[gl_SubgroupID] = max_value;
    }
	barrier();
	
	if (gl_SubgroupID == 0) {
    	max_value = data_max[gl_SubgroupInvocationID];
		subgroupBarrier();
        max_value = subgroupMax(max_value);
    }
	if (gl_LocalInvocationID.x == 0) {
		if(UV == 0){
			res[2 * gl_WorkGroupID.x + 0] = max_value;
		} else if(UV == 1){
			res[2 * gl_WorkGroupID.x + 1] = max_value;
		}
    }
	if(UV == 1 && idx == 0){
		cnt = cnt + 1;
	}
}