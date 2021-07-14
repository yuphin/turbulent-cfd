#version 450
#extension GL_GOOGLE_include_directive: require
#include "discretization.h"
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;
layout(constant_id = 0) const uint TURB_MODEL = 0;

layout(binding = 0) readonly buffer UVec
{
 	float U[];
};

layout(binding = 1) readonly buffer VVec
{
	float V[];
};

layout(binding = 21) buffer DT
{
	float dt;
};

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
layout(set = 1, binding = 17) buffer TURB_NU_I
{
    float NU_I[];
};
layout(set = 1, binding = 18) buffer TURB_NU_J
{
    float NU_J[];
};

layout(set = 1, binding = 19) buffer KOLD
{
    float K_old[];
};

layout(set = 1, binding = 20) buffer EPSOLD
{
    float EPS_old[];
};

float calculate_f1_sst(float omega, float dk_di, float dw_di, float k, float dist, float nu) {
    float cd_kw = max(2 * 0.856 * 1 / omega * dk_di * dw_di, 1e-10);
    float term = min(max(sqrt(k) / (0.09 * omega * dist), 500 * nu / (dist * dist * omega)),
                                    4 * 0.856 * k / (cd_kw * dist * dist));
    float pow_term = term * term;
    float q_term = pow_term * pow_term;
    float f1 = tanh(q_term);
    return f1;
}

float calculate_sst_term(float K_stencil[3], float EPS_stencil[3], float kij, float eij, float dist) {

    float dk_dx = (K_stencil[0] - K_stencil[1]) * inv_dx;
    float dw_dx = (EPS_stencil[0] - EPS_stencil[1]) * inv_dx;
    float dk_dy = (K_stencil[0] - K_stencil[2]) * inv_dy;
    float dw_dy = (EPS_stencil[0] - EPS_stencil[2]) * inv_dy;
    float f1_x = calculate_f1_sst(eij, dk_dx, dw_dx, kij, dist, nu);
    float f1_y = calculate_f1_sst(eij, dk_dy, dw_dy, kij, dist, nu);
    float res_x = 2 * (1 - f1_x) * 0.856 * 1 / eij * dk_dx * dw_dx;
    float res_y = 2 * (1 - f1_y) * 0.856 * 1 / eij * dk_dy * dw_dy;
    return res_x + res_y;
}


float convection_uKEPS(float U_stencil[2], float T_stencil[5]) {
    const int METHOD = 0;

    if (METHOD == 0) {
        float result = 0.0;
    
        float interp1 = (T_stencil[1] + T_stencil[0]) / 2;
        float interp2 = (T_stencil[2] + T_stencil[1]) / 2;
        float interp3 = (T_stencil[1] - T_stencil[0]) / 2;
        float interp4 = (T_stencil[2] - T_stencil[1]) /2;

        return inv_dx * ((U_stencil[1] * interp1 - U_stencil[0] * interp2) + 
                (abs(U_stencil[1]) * interp3 - abs(U_stencil[0]) * interp4));
    } else {
        // TODO
    }
}

float convection_vKEPS(float V_stencil[2], float T_stencil[5]) {
    const int METHOD = 0;
    if (METHOD == 0) {
        float interp1 = (T_stencil[1] + T_stencil[3]) / 2;
        float interp2 = (T_stencil[4] + T_stencil[1]) / 2;
        float interp3 = (T_stencil[1] - T_stencil[3]) / 2;
        float interp4 = (T_stencil[4] - T_stencil[1]) /2;

        return inv_dy * ((V_stencil[1] * interp1 - V_stencil[0] * interp2) + 
                (abs(V_stencil[1]) * interp3 - abs(V_stencil[0]) * interp4));
    } else {
        // TODO
    }
}

float laplacian_nu(float ar[5], float nu_i_stencil[2], float nu_j_stencil[2], float coeff) {
    // nu_i_stencilj[0] -> nu_i_stencilj(i,j)
    // nu_i_stencilj[1] -> nu_i_stencil(i-1,j) or nu_j_stencil(i, j-1)
    float inv_dx2 = inv_dx * inv_dx;
    float inv_dy2 = inv_dy * inv_dy;
    float i_diff = (nu + nu_i_stencil[0]) * (ar[0] - ar[1]) - (nu + nu_i_stencil[1]) * (ar[1] - ar[2]);
    float j_diff = (nu + nu_j_stencil[0]) * (ar[3] - ar[1]) - (nu + nu_j_stencil[1]) * (ar[1] - ar[4]);
    float result = 1 / coeff * (i_diff * inv_dx2 + j_diff * inv_dy2);
    return result;
}

float calculate_strain_tensor(float U_stencil[6], float V_stencil[6], float inv_dx, float inv_dy) {
    float shear_1 = (U_stencil[2] + U_stencil[3] - U_stencil[4] - U_stencil[5]) * (0.25 * inv_dy);
    float shear_2 = (V_stencil[2] + V_stencil[3] - V_stencil[4] - V_stencil[5]) * (0.25 * inv_dx);
    float shear = shear_1 + shear_2;
    return shear;
}

float mean_strain_rate_squared(float U_stencil[6], float V_stencil[6]) {
    // U_stencil offsets:
    // 0,0 / -1,0 / 0,1 / -1,1 / 0,-1 / -1,-1
    // V_stencil offsets:
    // 0,0 / 0,-1 / 1,0 / 1, -1 / -1,0 / -1,-1
    const int METHOD = 0;
    float result = 0;
    if (METHOD == 0) {
        float invdx2 = inv_dx * inv_dx;
        float invdy2 = inv_dy * inv_dy;
        float u_diff = U_stencil[0] - U_stencil[1];
        float v_diff = V_stencil[0] - V_stencil[1];
        float shear = calculate_strain_tensor(U_stencil, V_stencil, inv_dx, inv_dy);
        result = (u_diff * u_diff) * invdx2 + v_diff * v_diff * invdy2 + shear * shear;
    } else {
        // TODO
    }
    return result;
}

float mean_strain_rate_squared_out_shear(float U_stencil[6], float V_stencil[6], uint i, uint j, out float shear_out) {
    const int METHOD = 0;
    float result = 0;
    if (METHOD == 0) {
        float invdx2 = inv_dx * inv_dx;
        float invdy2 = inv_dy * inv_dy;
        float u_diff = U_stencil[0] - U_stencil[1];
        float v_diff = V_stencil[0] - V_stencil[1];
        float shear = calculate_strain_tensor(U_stencil, V_stencil, inv_dx, inv_dy);
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
    float is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 0){
        return;
    }
    float f2_coeff = 1;
    float nut = at(NU_T, i, j);
    float kij = at(K_old, i, j);
    float eij = at(EPS_old, i, j);
    float K_stencil[5] = float[5](at(K_old, i + 1, j), at(K_old, i, j), at(K_old, i - 1, j), at(K_old, i, j + 1),
                         at(K_old, i, j - 1));
    float EPS_stencil[5] = float[5](at(EPS_old, i + 1, j), at(EPS_old, i, j), at(EPS_old, i - 1, j), at(EPS_old, i, j + 1),
                           at(EPS_old, i, j - 1));
    float U_diff[2] = float[2](at(U, i - 1, j), at(U, i, j));
    float V_diff[2] = float[2](at(V, i, j - 1), at(V, i, j));
    float NU_I_diff[2] = float[2](at(NU_I, i, j), at(NU_I, i - 1, j));
    float NU_J_diff[2] = float[2](at(NU_J, i, j), at(NU_J, i, j - 1));
    float U_stencil[6] = float[6](at(U, i, j),         at(U, i - 1, j), at(U, i, j + 1),
                         at(U, i - 1, j + 1), at(U, i, j - 1), at(U, i - 1, j - 1));
    float V_stencil[6] = float[6](at(V, i, j),         at(V, i, j - 1), at(V, i + 1, j),
                         at(V, i + 1, j - 1), at(V, i - 1, j), at(V, i - 1, j - 1));
    float k1_1 = convection_uKEPS(U_diff, K_stencil);
    float k1_2 = convection_vKEPS(V_diff, K_stencil);
    float e1_1 = convection_uKEPS(U_diff, EPS_stencil);
    float e1_2 = convection_vKEPS(V_diff, EPS_stencil);
    float k2 = laplacian_nu(K_stencil, NU_I_diff, NU_J_diff, 1);
    float e2 = laplacian_nu(EPS_stencil, NU_I_diff, NU_J_diff, TURB_MODEL == 1 ? 1.3 : 1);
    float k3;
    if (TURB_MODEL != 3) {
        k3 = nut * mean_strain_rate_squared(U_stencil, V_stencil);
    } else {
        float shear_val = 0;
        k3 = nut * mean_strain_rate_squared_out_shear(U_stencil, V_stencil, i, j, shear_val);
        k3 = min(k3, 10 * 0.09 * kij * eij);
        // This is needed here for some reason
        if(k3 < 0) k3 = 0;
        at(S, i, j) = shear_val;
    }
    float e3 = (TURB_MODEL == 1 ? 1.44 : 5.0 / 9.0) * eij * k3 / kij;
    float e4 = TURB_MODEL == 1 ? f2_coeff * 1.92 * eij * eij / kij : 3.0 / 40 * eij * eij;
    float eij_mul = TURB_MODEL == 1 ? 1 : 0.09 * kij;
    float kij_new = kij + dt * (-(k1_1 + k1_2) + k2 + k3 - eij_mul * eij);
    float sst_term = 0;
    if (TURB_MODEL == 3) {
        float K_diff[3] = float[3](at(K_old, i, j), at(K_old, i - 1, j), at(K_old, i, j - 1));
        float EPS_diff[3] = float[3](at(EPS_old, i, j), at(EPS_old, i - 1, j), at(EPS_old, i, j - 1));
        float dist = at(dists, i, j);
        sst_term = calculate_sst_term(K_diff, EPS_diff, kij, eij, dist);
    }
    float epsij_new = eij + dt * (-(e1_1 + e1_2) + e2 + e3 - e4 + sst_term);
    at(K, i, j) = kij_new;
    at(EPS, i, j) = epsij_new;
}