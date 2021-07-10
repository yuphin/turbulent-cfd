#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
layout(local_size_x = 1024, local_size_y =1, local_size_z = 1) in;
layout(constant_id = 0) const uint KEPS = 0;
layout(binding = 0) buffer U
{
	double u[];
};

layout(binding = 1) buffer V
{
	double v[];
};

layout(binding = 14) buffer Residual //residual
{
	double res[];
};

layout(set = 1, binding = 13) buffer TURB_K
{
    double K[];
};
layout(set = 1, binding = 14) buffer TURB_EPS
{
    double EPS[];
};

layout(set = 1, binding = 22) buffer KepsResidual
{
	double kepsres[];
};

shared double data_maxu[32];
shared double data_maxv[32];
void main(){
	uint idx = gl_GlobalInvocationID.x;
	double max_value_u = 0;
	double max_value_v = 0;
	if(idx < size){
		if(KEPS == 0){
			max_value_u = abs(u[idx]);
			max_value_v = abs(v[idx]);
		} else if(KEPS == 1){
			max_value_u = K[idx];
			max_value_v = EPS[idx];
		}
	}

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
		if(KEPS == 0){
 			res[2*gl_WorkGroupID.x + 0] = max_value_u;
			res[2*gl_WorkGroupID.x + 1] = max_value_v;
		} else if(KEPS == 1){
			kepsres[2*gl_WorkGroupID.x + 0] = max_value_u;
			kepsres[2*gl_WorkGroupID.x + 1] = max_value_v;
		}
       
    }
}