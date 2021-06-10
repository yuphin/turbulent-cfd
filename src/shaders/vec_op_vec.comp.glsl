#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
layout(local_size_x = 1024, local_size_y = 1, local_size_z = 1) in;
const uint size_x = 102;
const uint size_y = 22;
const uint size = size_x * size_y;
layout(binding = 17) buffer V1
{
	float v1[];
};

layout(binding = 18) buffer V2
{
	float v2[];
};

layout(binding = 19) buffer Result
{
	float res[];
};

layout(constant_id = 0) const uint OP = 0;

void main() {
	uint idx = gl_GlobalInvocationID.x;
	// res = 0
	float sum = 0;
	if(idx < size){
		if(OP == 0) {
			res[idx] =  v1[idx] + v2[idx];
		
		}else if(OP == 1){
			res[idx] =  v1[idx] - v2[idx];
		}
	} 
}