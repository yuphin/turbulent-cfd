#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
layout(local_size_x = 1024, local_size_y = 1, local_size_z = 1) in;
const uint size_x = 102;
const uint size_y = 22;
const uint size = size_x * size_y;

layout(constant_id = 0) const uint MODE = 0;
layout(binding = 12) buffer V1 //d
{
	float v1[];
};

layout(binding = 13) buffer V2 //spmv
{
	float v2[];
};

layout(binding = 16) buffer Residual //residual
{
	float res[];
};

layout(binding = 19) buffer R // r
{
	float r[];
};

shared float data[32];

void main() {
	uint idx = gl_GlobalInvocationID.x;
	// res = 0
	float sum = 0;
	if(idx < size){
		if(MODE == 0){
			sum =  v1[idx] * v2[idx];
		} else if(MODE == 1){
			sum =  r[idx] * r[idx];
		}
		
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
	if (gl_LocalInvocationID.x == 0) {
		
        res[gl_WorkGroupID.x] = sum;
    }
}