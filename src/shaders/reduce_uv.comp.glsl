#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"

layout(local_size_x = 1024, local_size_y =1, local_size_z = 1) in;

layout(binding = 14) buffer Residual //residual
{
	float res[];
};

shared float data_maxu[32];
shared float data_maxv[32];
void main(){
	uint idx = gl_GlobalInvocationID.x;
	float max_value_u = 0;
	float max_value_v = 0;
	if(idx < num_wgs){
		max_value_u = res[2*idx + 0];
		max_value_v = res[2*idx + 1];
	}
	memoryBarrier();
	barrier();
	
	max_value_u = subgroupMax(max_value_u);
	max_value_v = subgroupMax(max_value_v);

	if (gl_SubgroupInvocationID == 0) {
        data_maxu[gl_SubgroupID] = max_value_u;
		data_maxv[gl_SubgroupID] = max_value_v;
    }
	barrier();
	
	if (gl_SubgroupID == 0) {
    	max_value_u = data_maxu[gl_SubgroupInvocationID];
		max_value_v = data_maxv[gl_SubgroupInvocationID];
		subgroupBarrier();
        max_value_u = subgroupMax(max_value_u);
		max_value_v = subgroupMax(max_value_v);
    }
	if (gl_LocalInvocationID.x == 0) {
		res[2 * gl_WorkGroupID.x + 0] = max_value_u;
		res[2 * gl_WorkGroupID.x + 1] = max_value_v;
    }
}