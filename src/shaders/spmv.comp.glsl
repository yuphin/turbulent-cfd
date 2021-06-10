#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
layout(local_size_x = 1024, local_size_y = 1, local_size_z = 1) in;
const uint size_x = 102;
const uint size_y = 22;
const uint num_rows = size_x * size_y;
const uint num_cols = size_x * size_y;
const uint num_diags = 5;

layout(binding = 10) buffer SparseMatrix
{
	float A[];
};


layout(binding = 11) buffer Offsets
{
	uint offsets[];
};

layout(binding = 12) buffer D
{
	float v[];
};

layout(binding = 13) buffer Result
{
	float res[];
};

void main() {
	uint row = gl_GlobalInvocationID.x;

	if(row < num_rows) {
		res[row] = 0;
		float sum = 0;
		for(int n = 0; n < num_diags; n++){
			uint col = row + offsets[n];

			float val = A[num_rows * n + row];

			if(col >= 0 && col < num_cols){
				sum += val * v[col];
			}	
		}
		res[row] += sum;
	}
}
