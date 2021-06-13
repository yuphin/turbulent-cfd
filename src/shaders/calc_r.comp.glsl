#version 450
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

layout(binding = 21) buffer DT
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

shared float data[32];

float laplacian(float ar[5]) {
    float inv_dx2 = inv_dx * inv_dx;
    float inv_dy2 = inv_dy * inv_dy;
    float result = (ar[0] - 2. * ar[1] + ar[2]) * inv_dx2 +
                   (ar[3] - 2. * ar[1] + ar[4]) * inv_dy2 ;
    return result;
}

void main() {
	uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
	float sum = 0.0;
    float is_fluid = at(cell_type, i, j);
	if((i < imax && j < jmax) && is_fluid == 1){
		float p_stencil[5] = float[5](at(p, i+1, j), at(p, i, j), 
								  at(p, i-1, j), at(p, i, j+1), 
								  at(p, i, j-1));
		float val = laplacian(p_stencil) - at(rs, i, j);
		sum = val * val;
    }

	sum = subgroupAdd(sum);
	if (gl_SubgroupInvocationID == 0) {
        data[gl_SubgroupID] = sum;
	
    }
	barrier();
	
    if (gl_SubgroupID == 0) {
        sum = data[gl_SubgroupInvocationID];
		subgroupBarrier();
        sum = subgroupAdd(sum); 
    }
	if (gl_SubgroupInvocationID == 0) {
        r[gl_WorkGroupID.x] = sum;
    }
  
}