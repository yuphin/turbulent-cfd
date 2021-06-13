#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
const uint size_x = 102;
const uint size_y = 22;
const uint size = size_x * size_y;
const uint dim = 32;

layout(local_size_x = 32, local_size_y =1, local_size_z = 1) in;

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


shared float data_minu[32];
shared float data_minv[32];
shared float data_maxu[32];
shared float data_maxv[32];

void main(){
	uint idx = gl_GlobalInvocationID.x;
	float max_value_u = 0;
	float min_value_u = 0;
	float max_value_v = 0;
	float min_value_v = 0;
	if(idx < dim){
		max_value_u = u_max[idx];
		min_value_u = u_min[idx];
		max_value_v = v_max[idx];
		min_value_v = v_min[idx];
	}

	max_value_u = subgroupMax(max_value_u);
	min_value_u = subgroupMin(min_value_u);
	max_value_v = subgroupMax(max_value_v);
	min_value_v = subgroupMin(min_value_v);

	// if (gl_SubgroupInvocationID == 0) {
    //    data_maxu[gl_SubgroupID] = max_value_u;
	// 	data_minu[gl_SubgroupID] = min_value_u;
	// 	data_maxv[gl_SubgroupID] = max_value_v;
	// 	data_minv[gl_SubgroupID] = min_value_v;
    // }
	// barrier();
	// if (gl_SubgroupID == 0) {
    // 	max_value_u = data_maxu[gl_SubgroupInvocationID];
	// 	min_value_u = data_minu[gl_SubgroupInvocationID];
	// 	max_value_v = data_maxv[gl_SubgroupInvocationID];
	// 	min_value_v = data_minv[gl_SubgroupInvocationID];
	// 	subgroupBarrier();
    //     max_value_u = subgroupMax(max_value_u);
	// 	min_value_u = subgroupMin(min_value_u);
	// 	max_value_v = subgroupMax(max_value_v);
	// 	min_value_v = subgroupMin(min_value_v);
    // }

	if (gl_LocalInvocationID.x == 0) {
		
        u_max[0] = max_value_u;
		u_min[0] = min_value_u;
		v_max[0] = max_value_v;
		v_min[0] = min_value_v;
    }
}