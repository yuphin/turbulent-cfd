#version 450
// Compile with glslangvalidator -V simulation.comp.glsl -o simulation.comp.spv
const  int imax = 102;
const  int jmax = 22;
const float nu = 0.1;
const float dx = 0.1;
const float dy = 0.1;
const float inv_dx = 10;
const float inv_dy = 10;
const float gamma = 0.5;
#define at(AR, I, J) AR[imax * (J) + (I)]
#define interpolate(A, i, j, i_offset, j_offset) (at(A, i, j) - at(A, i + i_offset, j + j_offset))/ 2
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;

layout(std430, binding = 0) readonly buffer U
{
 	float u[imax * jmax];
};

layout(std430, binding = 1) readonly buffer V
{
	float v[imax * jmax];
};

layout(std430, binding = 2) writeonly buffer F
{
	float f[imax * jmax];
};

layout(std430, binding = 3) writeonly buffer G
{
	float g[imax * jmax];
};

layout(binding = 4) uniform UBO
{
	float dt;
};

layout(binding = 5) readonly buffer CellType
{
	float cell_type[imax * jmax];
};


float convection_u(float U[6], float V[6]) {
    float result = 0.0;
    float interp1 = (U[0] + U[1]) / 2;
    float interp2 = (U[2] + U[1]) / 2;
    float interp3 = interp1 * 2;
    float interp4 = (U[1] - U[0]) / 2;
    float interp5 = interp2 * 2;
    float interp6 = (U[2] - U[1]) /2;

    float interp7 = (V[0] + V[1]) / 2;
    float interp8 = (U[3] + U[1]) / 2;
    float interp9 = (V[5] + V[4]) / 2;
    float interp10 = (U[4] + U[1]) /2;
    float interp11 = interp7 * 2;
    float interp12 = (U[1] - U[3]) /2;
    float interp13 = interp9 * 2;
    float interp14 = (U[4] - U[1]) /2;
   
    // dU^2/dx
    float result_fd = (interp1 * interp1 - interp2 * interp2) * inv_dx;
    float result_dc = gamma * (abs(interp3) / 2 * interp4 - abs(interp5) /2 * interp6) * inv_dx; 
    result += result_fd + result_dc;

    // dUV/dy
    result_fd = (interp7 * interp8 - interp9 * interp10) * inv_dy;
    result_dc = gamma * (abs(interp11) / 2 * interp12 -
                        abs(interp13) / 2 * interp14 ) * inv_dy;
    result += result_fd + result_dc;

    return result;

}
float convection_v(float U[6], float V[6]) {
    float result = 0.0;
    float interp1 = (V[3] + V[1]) / 2;
    float interp2 = (V[4] + V[1]) / 2;
    float interp3 = interp1 * 2;
    float interp4 = (V[1] - V[3]) / 2;
    float interp5 = interp2 * 2;
    float interp6 = (V[4] - V[1]) /2;

    float interp7 = (U[3] + U[1]) / 2;
    float interp8 = (V[0] + V[1]) / 2;
    float interp9 = (U[5] + U[2]) / 2;
    float interp10 = (V[2] + V[1]) /2;
    float interp11 = interp7 * 2;
    float interp12 = (V[1] - V[0]) /2;
    float interp13 = interp9 * 2;
    float interp14 = (V[2] - V[1]) /2;
   
    // dU^2/dx
    float result_fd = (interp1 * interp1 - interp2 * interp2) * inv_dy;
    float result_dc = gamma * (abs(interp3) / 2 * interp4 - abs(interp5) * interp6) * inv_dy; 
    result += result_fd + result_dc;

    // dUV/dy
    result_fd = (interp7 * interp8 - interp9 * interp10) * inv_dx;
    result_dc = gamma * (abs(interp11) / 2 * interp12 -
                        abs(interp13) / 2 * interp14 ) * inv_dx;
    result += result_fd + result_dc;
    return result;
}

float laplacian(float ar[6]) {
    float inv_dx2 = inv_dx * inv_dx;
    float inv_dy2 = inv_dy * inv_dy;
    float result = (ar[0] - 2. * ar[1] + ar[2]) * inv_dx2 +
                   (ar[3] - 2. * ar[1] + ar[4]) * inv_dy2 ;
    return result;
}

void main() {
    uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    float is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 0){
        return;
    }
    // 5-point + 1 stencil for U and V
    float u_stencil[6] = float[6](at(u, i+1, j), at(u, i, j), at(u, i-1, j), at(u, i, j+1), at(u, i, j-1), at(u, i-1, j+1));
    float v_stencil[6] = float[6](at(v, i+1, j), at(v, i, j), at(v, i-1, j), at(v, i, j+1), at(v, i, j-1), at(v, i+1, j-1));

   // Calculate fluxes
   at(f, i, j) = at(u, i, j ) + dt * (nu * laplacian(u_stencil) - convection_u(u_stencil,v_stencil));
   at(g, i, j) = at(v, i, j ) + dt * (nu * laplacian(v_stencil) - convection_v(u_stencil,v_stencil));
}