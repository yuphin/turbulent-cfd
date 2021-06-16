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
layout(binding = 21) buffer DT
{
	float dt;
};

void main(){
	float max_u = res[0];
	float max_v = res[1];
	float min_cond = min(dx / max_u, dy / max_v);
	if(nu != 0.0){
		float cond_spatial =  1.0 / (2.0 * nu) * ((dx2 * dy2) / (dx2 + dy2));
		min_cond = min(min_cond, cond_spatial);
	}

	dt = tau * min_cond;

}

