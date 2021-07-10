#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
layout(local_size_x = 1024, local_size_y = 1, local_size_z = 1) in;
 uint num_rows = imax * jmax;
 uint num_cols = imax * jmax;
layout(constant_id = 0) const uint MODE = 0;
layout(binding = 10) buffer SparseMatrix
{
	double A[];
};

layout(binding = 11) buffer Offsets
{
	uint a_offsets[];
};

layout(binding = 12) buffer D
{
	double v[];
};

layout(binding = 13) buffer Q
{
	double q[];
};

layout(binding = 15) buffer R
{
	double r[];
};

layout(binding = 17) buffer Z // z
{
	double z[];
};

layout(binding = 18) buffer MData // Preconditioner for A
{
	double M[];
};

layout(binding = 19) buffer MOffsets
{
	uint m_offsets[];
};

void main() {
	uint row = gl_GlobalInvocationID.x;

	if(row < num_rows) {
		if(MODE == 0){
			q[row] = 0;
			double sum = 0;
			for(int n = 0; n < 5; n++){
				uint col = row + a_offsets[n];
				double val = A[num_rows * n + row];
				if(col >= 0 && col < num_cols){
					sum += val * v[col];
				}	
			}
			q[row] += sum;
		} else if(MODE == 1){
			z[row] = 0;
			double sum = 0;
			for(int n = 0; n < num_diags; n++){
				uint col = row + m_offsets[n];
				double val = M[num_rows * n + row];
				if(col >= 0 && col < num_cols){
					sum += val * r[col];
				}	
			}
			z[row] += sum;
		}
	}
}
