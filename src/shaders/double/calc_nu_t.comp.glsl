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
    double kij = at(K, i, j);
    double epsij = at(EPS, i, j);
    if (TURB_MODEL == 1) {
        at(NU_T, i, j) = 0.09 * kij * kij / epsij + nu;
    } else if (TURB_MODEL == 2) {
        at(NU_T, i, j) = kij / epsij + nu;
    } else if (TURB_MODEL == 3) {
        const double a1 = 5.0 / 9.0;
        double dist = at(dists, i, j);
        double f2 = calculate_f2_sst(epsij, kij, dist, nu);
        at(NU_T, i, j) = a1 * kij / (max(a1 * epsij, at(S, i, j) * f2)) + nu;
    }
}