#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
const int dim = 32;

layout(local_size_x = 32, local_size_y =1, local_size_z = 1) in;
layout(constant_id = 0) const uint STORE = 0;
layout(binding = 14) buffer Buffer
{
	float r[];
};

layout(binding = 31) buffer Counter
{
	int cnt;
};

shared float data[32];

void main(){
	uint idx = gl_GlobalInvocationID.x;
	float sum = 0;
	if(idx < dim){
		
		sum = r[idx];
	}

	memoryBarrier();
	barrier();

	sum = subgroupAdd(sum);
	if (gl_LocalInvocationID.x == 0) {
		if(STORE == 0){
			r[gl_WorkGroupID.x] = sum;
		} else {
			int counter = cnt % 2;
			r[counter] = sum;
		}
	}
}
