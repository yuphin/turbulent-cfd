#version 450
#extension GL_GOOGLE_include_directive: require
#extension GL_KHR_shader_subgroup_vote: enable
const  int imax = 102;
const  int jmax = 22;
const float nu = 0.1;
const float dx = 0.1;
const float dy = 0.1;
const float inv_dx = 10;
const float inv_dy = 10;
const float inv_dx2 = inv_dx * inv_dx;
const float inv_dy2 = inv_dy * inv_dy;
const float gamma = 0.5;
const float omega = 1.0;
#define at(AR, I, J) AR[imax * (J) + (I)]
#define interpolate(A, i, j, i_offset, j_offset) (at(A, i, j) - at(A, i + i_offset, j + j_offset))/ 2
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;

layout(std430, binding = 0) buffer U
{
 	float u[imax * jmax];
};

layout(std430, binding = 1) buffer V
{
	float v[imax * jmax];
};

layout(std430, binding = 2) buffer F
{
	float f[imax * jmax];
};

layout(std430, binding = 3) buffer G
{
	float g[imax * jmax];
};

layout(binding = 4) uniform UBO
{
	float dt;
	float parity;
};

layout(binding = 5) readonly buffer CellType
{
	float cell_type[imax * jmax];
};

layout(binding = 6) buffer RS
{
	float rs[imax * jmax];
};

layout(binding = 7) buffer P
{
	float p[imax * jmax];
};

layout(constant_id = 0) const uint PARITY = 0;

float sor_helper(float p[5]) {
	float res = (p[0] + p[2]) * inv_dx2 + (p[3] + p[4]) * inv_dy2;
	return res;
}

void main() {
	uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    float is_fluid = at(cell_type, i, j);
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
	float p_stencil[5] = float[5](at(p, i+1, j), at(p, i, j), 
								  at(p, i-1, j), at(p, i, j+1), 
								  at(p, i, j-1));
	float coeff = omega / (2* (inv_dx2 + inv_dy2));
	at(p, i, j) = (1 - omega) * at(p, i, j) + coeff * (sor_helper(p_stencil) - (at(rs, i, j)));
}