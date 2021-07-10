#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
layout(constant_id = 0) const uint CALC_TEMP = 0;
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

	if(CALC_TEMP == 1){
		float inv_dx2 = inv_dx * inv_dx;
		float inv_dy2 = inv_dy * inv_dy;
		float cond_temp = 1 / (2* alpha * (inv_dx2 + inv_dy2));
		min_cond = min(min_cond, cond_temp);
	}

	dt = tau * min_cond;

}

