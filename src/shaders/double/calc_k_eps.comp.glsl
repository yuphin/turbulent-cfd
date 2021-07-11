#version 450
#extension GL_GOOGLE_include_directive: require
#include "discretization.h"
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;
layout(constant_id = 0) const uint TURB_MODEL = 0;

layout(binding = 0) readonly buffer UVec
{
 	double U[];
};

layout(binding = 1) readonly buffer VVec
{
	double V[];
};

layout(binding = 21) buffer DT
{
	double dt;
};

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

layout(set = 1, binding = 19) buffer KOLD
{
    double K_old[];
};

layout(set = 1, binding = 20) buffer EPSOLD
{
    double EPS_old[];
};

double calculate_f1_sst(double omega, double dk_di, double dw_di, double k, double dist, double nu) {
    double cd_kw = max(2 * 0.856 * 1 / omega * dk_di * dw_di, 1e-10);
    double term = min(max(sqrt(k) / (0.09 * omega * dist), 500 * nu / (dist * dist * omega)),
                                    4 * 0.856 * k / (cd_kw * dist * dist));
    double pow_term = term * term;
    double q_term = pow_term * pow_term;
    double f1 = dtanh(q_term);
    return f1;
}

double calculate_sst_term(double K_stencil[3], double EPS_stencil[3], double kij, double eij, double dist) {

    double dk_dx = (K_stencil[0] - K_stencil[1]) * inv_dx;
    double dw_dx = (EPS_stencil[0] - EPS_stencil[1]) * inv_dx;
    double dk_dy = (K_stencil[0] - K_stencil[2]) * inv_dy;
    double dw_dy = (EPS_stencil[0] - EPS_stencil[2]) * inv_dy;
    double f1_x = calculate_f1_sst(eij, dk_dx, dw_dx, kij, dist, nu);
    double f1_y = calculate_f1_sst(eij, dk_dy, dw_dy, kij, dist, nu);
    double res_x = 2 * (1 - f1_x) * 0.856 * 1 / eij * dk_dx * dw_dx;
    double res_y = 2 * (1 - f1_y) * 0.856 * 1 / eij * dk_dy * dw_dy;
    return res_x + res_y;
}


double convection_uKEPS(double U_stencil[2], double T_stencil[5]) {
    const int METHOD = 0;

    if (METHOD == 0) {
        double result = 0.0;
    
        double interp1 = (T_stencil[1] + T_stencil[0]) / 2;
        double interp2 = (T_stencil[2] + T_stencil[1]) / 2;
        double interp3 = (T_stencil[1] - T_stencil[0]) / 2;
        double interp4 = (T_stencil[2] - T_stencil[1]) /2;

        return inv_dx * ((U_stencil[1] * interp1 - U_stencil[0] * interp2) + 
                (abs(U_stencil[1]) * interp3 - abs(U_stencil[0]) * interp4));
    } else {
        // TODO
    }
}

double convection_vKEPS(double V_stencil[2], double T_stencil[5]) {
    const int METHOD = 0;
    if (METHOD == 0) {
        double interp1 = (T_stencil[1] + T_stencil[3]) / 2;
        double interp2 = (T_stencil[4] + T_stencil[1]) / 2;
        double interp3 = (T_stencil[1] - T_stencil[3]) / 2;
        double interp4 = (T_stencil[4] - T_stencil[1]) /2;

        return inv_dy * ((V_stencil[1] * interp1 - V_stencil[0] * interp2) + 
                gamma * (abs(V_stencil[1]) * interp3 - abs(V_stencil[0]) * interp4));
    } else {
        // TODO
    }
}

double laplacian_nu(double ar[5], double nu_i_stencil[2], double nu_j_stencil[2], double coeff) {
    // nu_i_stencilj[0] -> nu_i_stencilj(i,j)
    // nu_i_stencilj[1] -> nu_i_stencil(i-1,j) or nu_j_stencil(i, j-1)
    double inv_dx2 = inv_dx * inv_dx;
    double inv_dy2 = inv_dy * inv_dy;
    double i_diff = (nu + nu_i_stencil[0]) * (ar[0] - ar[1]) - (nu + nu_i_stencil[1]) * (ar[1] - ar[2]);
    double j_diff = (nu + nu_j_stencil[0]) * (ar[3] - ar[1]) - (nu + nu_j_stencil[1]) * (ar[1] - ar[4]);
    double result = 1 / coeff * (i_diff * inv_dx2 + j_diff * inv_dy2);
    return result;
}

double calculate_strain_tensor(double U_stencil[6], double V_stencil[6], double inv_dx, double inv_dy) {
    double shear_1 = (U_stencil[2] + U_stencil[3] - U_stencil[4] - U_stencil[5]) * (0.25 * inv_dy);
    double shear_2 = (V_stencil[2] + V_stencil[3] - V_stencil[4] - V_stencil[5]) * (0.25 * inv_dx);
    double shear = shear_1 + shear_2;
    return shear;
}

double mean_strain_rate_squared(double U_stencil[6], double V_stencil[6]) {
    // U_stencil offsets:
    // 0,0 / -1,0 / 0,1 / -1,1 / 0,-1 / -1,-1
    // V_stencil offsets:
    // 0,0 / 0,-1 / 1,0 / 1, -1 / -1,0 / -1,-1
    const int METHOD = 0;
    double result = 0;
    if (METHOD == 0) {
        double invdx2 = inv_dx * inv_dx;
        double invdy2 = inv_dy * inv_dy;
        double u_diff = U_stencil[0] - U_stencil[1];
        double v_diff = V_stencil[0] - V_stencil[1];
        double shear = calculate_strain_tensor(U_stencil, V_stencil, inv_dx, inv_dy);
        result = (u_diff * u_diff) * invdx2 + v_diff * v_diff * invdy2 + shear * shear;
    } else {
        // TODO
    }
    return result;
}

double mean_strain_rate_squared_out_shear(double U_stencil[6], double V_stencil[6], uint i, uint j, out double shear_out) {
    const int METHOD = 0;
    double result = 0;
    if (METHOD == 0) {
        double invdx2 = inv_dx * inv_dx;
        double invdy2 = inv_dy * inv_dy;
        double u_diff = U_stencil[0] - U_stencil[1];
        double v_diff = V_stencil[0] - V_stencil[1];
        double shear = calculate_strain_tensor(U_stencil, V_stencil, inv_dx, inv_dy);
        result = (u_diff * u_diff) * invdx2 + v_diff * v_diff * invdy2 + shear * shear;
        shear_out = shear;
    } else {
        // TODO
    }
    return result;
}

void main() {
    uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    double is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 0){
        return;
    }
    double f2_coeff = 1;
    double nut = at(NU_T, i, j);
    double kij = at(K_old, i, j);
    double eij = at(EPS_old, i, j);
    double K_stencil[5] = double[5](at(K_old, i + 1, j), at(K_old, i, j), at(K_old, i - 1, j), at(K_old, i, j + 1),
                         at(K_old, i, j - 1));
    double EPS_stencil[5] = double[5](at(EPS_old, i + 1, j), at(EPS_old, i, j), at(EPS_old, i - 1, j), at(EPS_old, i, j + 1),
                           at(EPS_old, i, j - 1));
    double U_diff[2] = double[2](at(U, i - 1, j), at(U, i, j));
    double V_diff[2] = double[2](at(V, i, j - 1), at(V, i, j));
    double NU_I_diff[2] = double[2](at(NU_I, i, j), at(NU_I, i - 1, j));
    double NU_J_diff[2] = double[2](at(NU_J, i, j), at(NU_J, i, j - 1));
    double U_stencil[6] = double[6](at(U, i, j),         at(U, i - 1, j), at(U, i, j + 1),
                         at(U, i - 1, j + 1), at(U, i, j - 1), at(U, i - 1, j - 1));
    double V_stencil[6] = double[6](at(V, i, j),         at(V, i, j - 1), at(V, i + 1, j),
                         at(V, i + 1, j - 1), at(V, i - 1, j), at(V, i - 1, j - 1));
    double k1_1 = convection_uKEPS(U_diff, K_stencil);
    double k1_2 = convection_vKEPS(V_diff, K_stencil);
    double e1_1 = convection_uKEPS(U_diff, EPS_stencil);
    double e1_2 = convection_vKEPS(V_diff, EPS_stencil);
    double k2 = laplacian_nu(K_stencil, NU_I_diff, NU_J_diff, 1);
    double e2 = laplacian_nu(EPS_stencil, NU_I_diff, NU_J_diff, TURB_MODEL == 1 ? 1.3 : 1);
    double k3;
    if (TURB_MODEL != 3) {
        k3 = nut * mean_strain_rate_squared(U_stencil, V_stencil);
    } else {
        double shear_val = 0;
        k3 = nut * mean_strain_rate_squared_out_shear(U_stencil, V_stencil, i, j, shear_val);
        k3 = min(k3, 10 * 0.09 * kij * eij);
        // This is needed here for some reason
        if(k3 < 0) k3 = 0;
        at(S, i, j) = shear_val;
    }
    double e3 = (TURB_MODEL == 1 ? 1.44 : 5.0 / 9.0) * eij * k3 / kij;
    double e4 = TURB_MODEL == 1 ? f2_coeff * 1.92 * eij * eij / kij : 3.0 / 40 * eij * eij;
    double eij_mul = TURB_MODEL == 1 ? 1 : 0.09 * kij;
    double kij_new = kij + dt * (-(k1_1 + k1_2) + k2 + k3 - eij_mul * eij);
    double sst_term = 0;
    if (TURB_MODEL == 3) {
        double K_diff[3] = double[3](at(K_old, i, j), at(K_old, i - 1, j), at(K_old, i, j - 1));
        double EPS_diff[3] = double[3](at(EPS_old, i, j), at(EPS_old, i - 1, j), at(EPS_old, i, j - 1));
        double dist = at(dists, i, j);
        sst_term = calculate_sst_term(K_diff, EPS_diff, kij, eij, dist);
    }
    double epsij_new = eij + dt * (-(e1_1 + e1_2) + e2 + e3 - e4 + sst_term);
    at(K, i, j) = kij_new;
    at(EPS, i, j) = epsij_new;
}