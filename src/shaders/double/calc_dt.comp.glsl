#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
layout(constant_id = 0) const uint CALC_TEMP = 0;
layout(constant_id = 1) const uint TURB_MODEL = 0;
layout(local_size_x = 1, local_size_y =1, local_size_z = 1) in;

layout(set = 1, binding = 21) buffer NutResidual
{
	double nures[];
};

layout(set = 1, binding = 22) buffer KepsResidual
{
	double kepsres[];
};
layout(binding = 14) buffer Residual
{
	double res[];
};
layout(binding = 21) buffer DT
{
	double dt;
};

void main(){
	double max_u = res[0];
	double max_v = res[1];
	double min_cond = min(dx / max_u, dy / max_v);
	double nu_min = nures[0];
	nu_min = nu_min == 0 ? nu : nu_min;
	if(nu_min != 0.0){
		double cond_spatial =  1.0 / (2.0 * nu) * ((dx2 * dy2) / (dx2 + dy2));
		min_cond = min(min_cond, cond_spatial);
	}
	
	if(CALC_TEMP == 1){
		double inv_dx2 = inv_dx * inv_dx;
		double inv_dy2 = inv_dy * inv_dy;
		double cond_temp = 1 / (2* alpha * (inv_dx2 + inv_dy2));
		min_cond = min(min_cond, cond_temp);
	}
	if(TURB_MODEL != 0){
		double k_max = kepsres[0];
		double eps_max = kepsres[1];
  		double cond_5 = 1 / (2 * k_max * (1 / dx2 + 1 / dy2));
        double cond_6;
        if (TURB_MODEL == 1) {
            cond_6 = 1 / (2 * eps_max * (1 / dx2 + 1 / dy2));
        } else if (TURB_MODEL == 2 || TURB_MODEL == 3) {
            cond_6 = 1 / (2 * (eps_max * 0.09 * k_max) * (1 / dx2 + 1 / dy2));
        }
        min_cond = min(min_cond, cond_5);
        min_cond = min(min_cond, cond_6);
	}

	dt = tau * min_cond;

}

