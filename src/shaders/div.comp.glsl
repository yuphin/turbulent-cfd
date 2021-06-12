#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable

layout(constant_id = 0) const uint STORE = 0;
layout(local_size_x = 1, local_size_y =1, local_size_z = 1) in;

layout(binding = 14) buffer R
{
	float r[];
};

layout(binding = 30) buffer Deltas
{
	float d[];
};

layout(binding = 31) buffer Counter
{
	int cnt;
};

void main(){
	int counter = cnt % 2;
	if(STORE == 0){
		r[counter] =  d[0] / r[counter];
	} else {
		d[1] = d[0]; // old
		d[0] = r[counter]; // new
		r[counter] = r[counter] / d[1];
	}

}

