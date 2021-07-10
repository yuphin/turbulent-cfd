#version 450
#extension GL_GOOGLE_include_directive: require
#include "discretization.h"
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;
layout(constant_id = 0) const uint TURB_MODEL = 0;
layout(binding = 5) readonly buffer CellType
{
	double cell_type[];
};

layout(set = 1, binding = 12) buffer TURB_NU_T
{
    double NU_T[];
};
layout(set = 1, binding = 13) buffer TURB_K
{
    double K[];
};
layout(set = 1, binding = 14) buffer TURB_EPS
{
    double EPS[];
};
layout(set = 1, binding = 15) buffer Distances
{
    double dists[];
};
layout(set = 1, binding = 16) buffer StrainTensor
{
    double S[];
};
layout(set = 1, binding = 17) buffer TURB_NU_I
{
    double NU_I[];
};
layout(set = 1, binding = 18) buffer TURB_NU_J
{
    double NU_J[];
};

double calculate_f2_sst(double omega, double k, double dist, double nu) {
    double max_sqr = max(2 * sqrt(k) / (0.09 * omega * dist), 500 * nu / (dist * dist * omega));
    return dtanh(max_sqr * max_sqr);
}

void main() {
    uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    double is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 0){
        return;
    }
    double num_i = (at(K, i, j) + at(K, i + 1, j)) / 2;
    double denom_i = (at(EPS, i, j) + at(EPS, i + 1, j)) / 2;

    double num_j = (at(K, i, j) + at(K, i, j + 1)) / 2;
    double denom_j = (at(EPS, i, j) + at(EPS, i, j + 1)) / 2;
    if (TURB_MODEL == 1) {
        at(NU_I, i, j) = 0.09 * num_i * num_i / denom_i;
        at(NU_J, i, j) = 0.09 * num_j * num_j / denom_j;
    } else if (TURB_MODEL == 2) {
        at(NU_I, i, j) = 0.5 * num_i / denom_i;
        at(NU_J, i, j) = 0.5 * num_j / denom_j;
    } else if (TURB_MODEL == 3) {
        double a1 = 5.0 / 9.0;
        double dist = at(dists, i, j);
        double f2_1 = calculate_f2_sst(denom_i, num_i, dist, nu);
        double f2_2 = calculate_f2_sst(denom_j, num_j, dist, nu);
        at(NU_I, i, j) = 0.85 * a1 * num_i / max(a1 * denom_i, at(S, i, j) * f2_1);
        at(NU_J, i, j) = 0.5 * a1 * num_j / max(a1 * denom_j, at(S, i, j) * f2_2);
    }
}