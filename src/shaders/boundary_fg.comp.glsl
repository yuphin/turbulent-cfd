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
const float UIN_2 = 1;
const float VIN_2 = 0;
const float UI = 0;
const float VI = 0;
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


layout(binding = 5) readonly buffer CellType
{
	float cell_type[imax * jmax];
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

	switch(neighbors){
			case 0: 
			{
				at(f, i, j) = at(u, i, j);
			}
			break;
			case 1: 
			{
				at(f, i-1, j) = at(u, i-1, j);
			}
			break;
			case 2: 
			{
				at(g, i, j) = at(v, i, j);
			}
			break;
			case 3: 
			{
				at(g, i, j-1) = at(v, i, j-1);
			}
			break;
	}
}