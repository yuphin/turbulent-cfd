#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable

layout(constant_id = 0) const uint STORE = 0;
layout(local_size_x = 1024, local_size_y = 1, local_size_z = 1) in;
const uint size_x = 102;
const uint size_y = 22;
const uint size = size_x * size_y;
layout(binding = 16) buffer P
{
	float p[];
};


void main(){
    uint idx = gl_GlobalInvocationID.x;
	if(idx < size){
        p[idx] = -p[idx];
    }
}

