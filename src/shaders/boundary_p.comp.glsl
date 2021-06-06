#version 450
#extension GL_GOOGLE_include_directive: require
#extension GL_KHR_shader_subgroup_arithmetic : enable
const  int imax = 102;
const  int jmax = 22;
const float nu = 0.1;
const float dx = 0.1;
const float dy = 0.1;
const float inv_dx = 10;
const float inv_dy = 10;
const float gamma = 0.5;
const float omega = 1.7;
const float PI = 0.0;
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

layout(binding = 8) buffer R
{
	float r[];
};

layout(binding = 9) buffer Neighborhood {
    uint neighborhood[];
};


void main() {
	uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    float is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 1) {
        return;
    }
	uint type = at(neighborhood, i, j) >> 8;
	uint neighbors = at(neighborhood, i, j) & 0xFF;
	if(neighbors == 0xFF) {
		return;
	}
	if(type == 0) {
		// Outlet
		at(p, i, j) = PI;
		//at(p, i, j) = imax * (j) + (i);
	} else  {
		// No-slip
		float p_stencil[8] = float[8](at(p, i+1, j), at(p, i - 1, j), 
								  at(p, i, j + 1), at(p, i, j-1), 
								  (at(p, i+1, j) + at(p, i, j + 1)) / 2,
								  (at(p, i+1, j) + at(p, i, j - 1)) / 2,
								  (at(p, i-1, j) + at(p, i, j + 1)) / 2,
								  (at(p, i-1, j) + at(p, i, j - 1)) / 2);
		at(p, i, j) = p_stencil[neighbors];
		//at(p, i, j) = imax * (j) + (i);//p_stencil[neighbors];
	}
}