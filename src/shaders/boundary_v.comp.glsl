#version 450
#extension GL_GOOGLE_include_directive: require
#extension GL_KHR_shader_subgroup_arithmetic : enable
const  int imax = 102;
const  int jmax = 22;
const float nu = 0.1;
const float dx = 0.1;
const float dy = 0.1;
const float inv_dx = 10;
const float inv_dy = 10;
const float gamma = 0.5;
const float omega = 1.7;
const float PI = 0.0;
const float UIN_2 = 1;
const float VIN_2 = 0;
const float UI = 0;
const float VI = 0;
#define at(AR, I, J) AR[imax * (J) + (I)]
#define interpolate(A, i, j, i_offset, j_offset) (at(A, i, j) - at(A, i + i_offset, j + j_offset))/ 2
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;

layout(std430, binding = 0) buffer U
{
 	float u[imax * jmax];
};

layout(std430, binding = 1) buffer V
{
	float v[imax * jmax];
};


layout(binding = 5) readonly buffer CellType
{
	float cell_type[imax * jmax];
};


layout(binding = 9) buffer Neighborhood {
    uint neighborhood[];
};


void main() {
	uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    float is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 1) {
        return;
    }
	uint type = at(neighborhood, i, j) >> 8;
	uint neighbors = at(neighborhood, i, j) & 0xFF;
	if(neighbors == 0xFF) {
		return;
	}

	switch(type) {
		case 0:
			// Outlet
			break;
		case 1:
			// Inlet
			switch(neighbors){
				case 0: 
				{
					at(u, i, j) = UIN_2;
					at(v, i, j) = 2 * VIN_2 - at(v, i+1, j);
				}
				break;
				case 1: 
				{
					at(u, i - 1, j) = UIN_2;
					at(v, i, j) = 2 * VIN_2 - at(v, i-1, j);
				}
				break;
				case 2: 
				{
					at(u, i, j) = 2 * UIN_2 - at(u, i, j+1);
					at(v, i, j) = VIN_2;
				}
				break;
				case 3: 
				{
					at(u, i, j) = 2 * UIN_2 - at(u, i, j-1);
					at(v, i, j-1) = VIN_2;
				
				}
				break;
			}
			break;
		case 2:
		{
			// NoSlip
			switch(neighbors){
				case 0: 
				{
					at(u, i, j) = 0;
					at(v, i, j) = 2 * VI - at(v, i+1, j);
				}
				break;
				case 1: 
				{
					at(u, i-1, j) = 0;
					at(v, i, j) = 2 * VI - at(v, i-1, j);
				}
				break;
				case 2: 
				{
					at(v, i, j) = 0;
					at(u, i, j) = 2 * UI - at(u, i, j+1);
				
				}
				break;
				case 3: 
				{
					at(v, i, j-1) = 0;
					at(u, i, j) = 2 * UI - at(u, i, j-1);
				
				}
				break;
				case 4: 
				{
					at(u, i, j) = 0;
					at(v, i, j) = 0;
					at(u, i-1, j) = 2 * UI - at(u, i-1, j+1);
					at(v, i, j-1) = 2 * VI - at(v, i+1, j-1);
				
				}
				break;
				case 5: 
				{
					at(u, i, j) = 0;
					at(v, i, j) = 2 * VI - at(v, i+1, j);
					at(u, i-1, j) = 2 * UI - at(u, i-1, j-1);
					at(v, i, j-1) = 0;
				
				}
				break;
				case 6: 
				{
					at(u, i, j) = 2 * UI - at(u, i, j+1);
					at(v, i, j) = 0;
					at(u, i-1, j) = 0;
					at(v, i, j-1) = 2 * VI - at(v, i-1, j-1);
				
				}
				break;
				case 7: 
				{
					at(u, i, j) = 2 * UI - at(u, i, j-1);
					at(v, i, j) = 2 * VI - at(v, i-1, j);
					at(u, i-1, j) = 0;
					at(v, i, j-1) = 0;
				}
				break;
			}
		}	
			break;
		case 3:
			// FreeSlip	
			switch(neighbors){
				case 0: 
				{
					at(u, i, j) = 0;
				}
				break;
				case 1: 
				{
					at(u, i-1, j) = 0;
				}
				break;
				case 2: 
				{
					at(v, i, j) = 0;
				
				}
				break;
				case 3: 
				{
					at(v, i, j-1) = 0;
				
				}
				break;
			}
			break;
		 default:
			break;

	}
}