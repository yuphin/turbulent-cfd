#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
const int dim = 102 * 22;

layout(local_size_x = 32, local_size_y =1, local_size_z = 1) in;

layout(binding = 30) buffer Counter
{
	int cnt;
};

layout(binding = 31) buffer Buffer
{
	float r[];
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

	if(gl_LocalInvocationID.x == 0){
		atomicAdd(cnt, 1);
	}
}

