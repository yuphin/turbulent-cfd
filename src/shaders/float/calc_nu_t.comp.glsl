#version 450
#extension GL_GOOGLE_include_directive: require
#include "discretization.h"
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;
layout(constant_id = 0) const uint TURB_MODEL = 0;
layout(binding = 5) readonly buffer CellType
{
	float cell_type[];
};

layout(set = 1, binding = 12) buffer TURB_NU_T
{
    float NU_T[];
};
layout(set = 1, binding = 13) buffer TURB_K
{
    float K[];
};
layout(set = 1, binding = 14) buffer TURB_EPS
{
    float EPS[];
};
layout(set = 1, binding = 15) buffer Distances
{
    float dists[];
};
layout(set = 1, binding = 16) buffer StrainTensor
{
    float S[];
};

float calculate_f2_sst(float omega, float k, float dist, float nu) {
    float max_sqr = max(2 * sqrt(k) / (0.09 * omega * dist), 500 * nu / (dist * dist * omega));
    return tanh(max_sqr * max_sqr);
}

void main() {
    uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    float is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 0){
        return;
    }
    float kij = at(K, i, j);
    float epsij = at(EPS, i, j);
    if (TURB_MODEL == 1) {
        at(NU_T, i, j) = 0.09 * kij * kij / epsij + nu;
    } else if (TURB_MODEL == 2) {
        at(NU_T, i, j) = kij / epsij + nu;
    } else if (TURB_MODEL == 3) {
        const float a1 = 5.0 / 9.0;
        float dist = at(dists, i, j);
        float f2 = calculate_f2_sst(epsij, kij, dist, nu);
        at(NU_T, i, j) = a1 * kij / (max(a1 * epsij, at(S, i, j) * f2)) + nu;
    }
}