#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable

layout(constant_id = 0) const uint STORE = 0;
layout(local_size_x = 1, local_size_y =1, local_size_z = 1) in;

layout(binding = 31) buffer Counter
{
	int cnt;
};

void main(){
	cnt = cnt + 1;
}

