#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
layout(local_size_x = 1024, local_size_y = 1, local_size_z = 1) in;
const uint size_x = 102;
const uint size_y = 22;
const uint size = size_x * size_y;
layout(binding = 20) uniform a
{
	float alpha;
};
layout(binding = 19) buffer V1
{
	float v1[];
};

layout(binding = 21) buffer V2
{
	float v2[];
};

layout(binding = 22) buffer Result
{
	float res[];
};

void main() {
	uint idx = gl_GlobalInvocationID.x;
	// res = 0
	if(idx < size){
		res[idx] = v1[idx] + alpha * v2[idx];
	} 
}