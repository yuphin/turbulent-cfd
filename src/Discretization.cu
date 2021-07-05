#pragma once
#ifndef __CUDACC__
#define __CUDACC__
#endif
#include "Utilities.hpp"
#include <cuda.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#define at(AR, I, J) AR[imax * (J) + (I)]

__device__ Real convection_u(Real U[6], Real V[6], Real inv_dx, Real inv_dy, Real gamma) {
    Real result = 0.0;
    Real interp1 = (U[0] + U[1]) / 2;
    Real interp2 = (U[2] + U[1]) / 2;
    Real interp3 = interp1 * 2;
    Real interp4 = (U[1] - U[0]) / 2;
    Real interp5 = interp2 * 2;
    Real interp6 = (U[2] - U[1]) / 2;

    Real interp7 = (V[0] + V[1]) / 2;
    Real interp8 = (U[3] + U[1]) / 2;
    Real interp9 = (V[5] + V[4]) / 2;
    Real interp10 = (U[4] + U[1]) / 2;
    Real interp11 = interp7 * 2;
    Real interp12 = (U[1] - U[3]) / 2;
    Real interp13 = interp9 * 2;
    Real interp14 = (U[4] - U[1]) / 2;

    // dU^2/dx
    Real result_fd = (interp1 * interp1 - interp2 * interp2) * inv_dx;
    Real result_dc = gamma * (abs(interp3) / 2 * interp4 - abs(interp5) / 2 * interp6) * inv_dx;
    result += result_fd + result_dc;

    // dUV/dy
    result_fd = (interp7 * interp8 - interp9 * interp10) * inv_dy;
    result_dc = gamma * (abs(interp11) / 2 * interp12 - abs(interp13) / 2 * interp14) * inv_dy;
    result += result_fd + result_dc;

    return result;
}

__device__  Real convection_uT(Real U[2], Real T[5], Real inv_dx, Real gamma) {
    Real result = 0.0;

    Real interp1 = (T[1] + T[0]) / 2;
    Real interp2 = (T[2] + T[1]) / 2;
    Real interp3 = (T[1] - T[0]) / 2;
    Real interp4 = (T[2] - T[1]) / 2;

    return inv_dx * ((U[1] * interp1 - U[0] * interp2) + gamma * (abs(U[1]) * interp3 - abs(U[0]) * interp4));
}

__device__  Real convection_vT(Real V[2], Real T[5], Real inv_dy, Real gamma) {
    Real result = 0.0;

    Real interp1 = (T[1] + T[3]) / 2;
    Real interp2 = (T[4] + T[1]) / 2;
    Real interp3 = (T[1] - T[3]) / 2;
    Real interp4 = (T[4] - T[1]) / 2;

    return inv_dy * ((V[1] * interp1 - V[0] * interp2) + gamma * (abs(V[1]) * interp3 - abs(V[0]) * interp4));
}

__device__  Real convection_v(Real U[6], Real V[6], Real inv_dx, Real inv_dy, Real gamma) {
    Real result = 0.0;
    Real interp1 = (V[3] + V[1]) / 2;
    Real interp2 = (V[4] + V[1]) / 2;
    Real interp3 = interp1 * 2;
    Real interp4 = (V[1] - V[3]) / 2;
    Real interp5 = interp2 * 2;
    Real interp6 = (V[4] - V[1]) / 2;

    Real interp7 = (U[3] + U[1]) / 2;
    Real interp8 = (V[0] + V[1]) / 2;
    Real interp9 = (U[5] + U[2]) / 2;
    Real interp10 = (V[2] + V[1]) / 2;
    Real interp11 = interp7 * 2;
    Real interp12 = (V[1] - V[0]) / 2;
    Real interp13 = interp9 * 2;
    Real interp14 = (V[2] - V[1]) / 2;

    // dU^2/dx
    Real result_fd = (interp1 * interp1 - interp2 * interp2) * inv_dy;
    Real result_dc = gamma * (abs(interp3) / 2 * interp4 - abs(interp5) * interp6) * inv_dy;
    result += result_fd + result_dc;

    // dUV/dy
    result_fd = (interp7 * interp8 - interp9 * interp10) * inv_dx;
    result_dc = gamma * (abs(interp11) / 2 * interp12 - abs(interp13) / 2 * interp14) * inv_dx;
    result += result_fd + result_dc;
    return result;
}

__device__  Real laplacian(Real ar[6], Real inv_dx, Real inv_dy) {
    Real inv_dx2 = inv_dx * inv_dx;
    Real inv_dy2 = inv_dy * inv_dy;
    Real result = (ar[0] - 2. * ar[1] + ar[2]) * inv_dx2 + (ar[3] - 2. * ar[1] + ar[4]) * inv_dy2;
    return result;
}

__device__  Real laplacian_5(Real ar[5], Real inv_dx, Real inv_dy) {
    Real inv_dx2 = inv_dx * inv_dx;
    Real inv_dy2 = inv_dy * inv_dy;
    Real result = (ar[0] - 2. * ar[1] + ar[2]) * inv_dx2 + (ar[3] - 2. * ar[1] + ar[4]) * inv_dy2;
    return result;
}