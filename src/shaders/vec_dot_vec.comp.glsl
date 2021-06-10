#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
layout(local_size_x = 1024, local_size_y = 1, local_size_z = 1) in;
const uint size_x = 102;
const uint size_y = 22;
const uint size = size_x * size_y;
layout(binding = 14) buffer V1
{
	float v1[];
};

layout(binding = 15) buffer V2
{
	float v2[];
};

layout(binding = 16) buffer Result
{
	float res[];
};

shared float data[32];

void main() {
	uint idx = gl_GlobalInvocationID.x;
	// res = 0
	float sum = 0;
	if(idx < size){
		//data[tid] = v1[idx] * v2[idx];
		sum =  v1[idx] * v2[idx];
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
        res[gl_WorkGroupID.x] = sum;
    }
	
//	for(uint i = gl_WorkGroupSize.x / 2 ; i >= 1; i = i/2){
//		barrier();
//		if(gl_LocalInvocationID.x < i) {
//			data[tid] += data[tid + i];
//		}
//	}

//	if(tid == 0) {
//		res[gl_WorkGroupID.x] = data[tid];
//	}



	
	
}