#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable

layout(constant_id = 0) const uint STORE = 0;
layout(local_size_x = 1, local_size_y =1, local_size_z = 1) in;
const float dx = 0.1;
const float dy = 0.1;
const float dx2 = dx * dx;
const float dy2 = dy * dy;
const float tau = 0.5;
const float nu = 0.1;
layout(binding = 17) buffer UMaxOut
{
	float u_max[];
};

layout(binding = 18) buffer UMinOut
{
	float u_min[];
};

layout(binding = 19) buffer VMaxOut
{
	float v_max[];
};

layout(binding = 20) buffer VMinOut
{
	float v_min[];
};

layout(binding = 21) buffer DT
{
	float dt;
};

void main(){
	float max_u = max(abs(u_min[0]), abs(u_max[0]));
	float max_v = max(abs(u_min[0]), abs(u_max[0]));
	float min_cond = (dx / max_u, dy / max_v);
	if(nu != 0.0){
		float cond_spatial =  1.0 / (2.0 * nu) * ((dx2 * dy2) / (dx2 + dy2));
		min_cond = min(min_cond, cond_spatial);
	}

	dt = tau * min_cond;

}

