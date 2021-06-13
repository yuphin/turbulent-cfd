#version 450
#extension GL_GOOGLE_include_directive: require

const  int imax = 102;
const  int jmax = 22;
const float nu = 0.1;
const float dx = 0.1;
const float dy = 0.1;
const float inv_dx = 10;
const float inv_dy = 10;
const float gamma = 0.5;
#define at(AR, I, J) AR[imax * (J) + (I)]
#define interpolate(A, i, j, i_offset, j_offset) (at(A, i, j) - at(A, i + i_offset, j + j_offset))/ 2
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;

layout(std430, binding = 0) readonly buffer U
{
 	float u[imax * jmax];
};

layout(std430, binding = 1) readonly buffer V
{
	float v[imax * jmax];
};

layout(std430, binding = 2) writeonly buffer F
{
	float f[imax * jmax];
};

layout(std430, binding = 3) writeonly buffer G
{
	float g[imax * jmax];
};

layout(binding = 21) buffer DT
{
	float dt;
};

layout(binding = 5) readonly buffer CellType
{
	float cell_type[imax * jmax];
};

layout(binding = 6) writeonly buffer RS
{
	float rs[imax * jmax];
};

float calculate_rs(float f[2], float g[2]) {
	float f_diff = inv_dx * ( f[0] - f[1]);
	float g_diff = inv_dy * ( g[0] - g[1]);
	return (f_diff + g_diff) * 1/ dt;
}


void main() {
    uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    float is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 0){
        return;
    }
	float f_diff[2] = float[2](at(f, i, j), at(f, i-1, j));
	float g_diff[2] = float[2](at(g, i, j), at(g, i, j - 1));
	at(rs, i, j) = calculate_rs(f_diff, g_diff);
	
}