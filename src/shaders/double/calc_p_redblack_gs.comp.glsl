#version 450
#extension GL_KHR_shader_subgroup_vote: enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
#define at(AR, I, J) AR[imax * (J) + (I)]
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;

layout(binding = 5) readonly buffer CellType
{
	double cell_type[];
};

layout(binding = 6) buffer RS
{
	double rs[];
};

layout(binding = 7) buffer P
{
	double p[];
};

layout(constant_id = 0) const uint PARITY = 0;

double sor_helper(double p[5], double inv_dx2, double inv_dy2) {

	double res = (p[0] + p[2]) * inv_dx2 + (p[3] + p[4]) * inv_dy2;
	return res;
}

void main() {
	uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    double is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 0) {
        return;
    }
	if(PARITY == 0) {
		// Red 
		if(((i + j) % 2) == 0){
			return;
		}
	} else if(PARITY == 1){
		// Black
		if(((i + j) % 2) == 1){
			return;
		}
	}
	const double omega = 1.7;
	double inv_dx2 = inv_dx * inv_dx;
	double inv_dy2 = inv_dy * inv_dy;
	double p_stencil[5] = double[5](at(p, i+1, j), at(p, i, j), 
								  at(p, i-1, j), at(p, i, j+1), 
								  at(p, i, j-1));
	double coeff = omega / (2* (inv_dx2 + inv_dy2));
	at(p, i, j) = (1- omega) * at(p,i,j) + 
			      coeff * (sor_helper(p_stencil, inv_dx2, inv_dy2) - (at(rs, i, j)));
}