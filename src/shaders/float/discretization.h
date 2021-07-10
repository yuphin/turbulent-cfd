#define at(AR, I, J) AR[imax * (J) + (I)]
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
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

float convection_uT(float U[2], float T[5]) {
    float result = 0.0;
    
    float interp1 = (T[1] + T[0]) / 2;
    float interp2 = (T[2] + T[1]) / 2;
    float interp3 = (T[1] - T[0]) / 2;
    float interp4 = (T[2] - T[1]) /2;

    return inv_dx * ((U[1] * interp1 - U[0] * interp2) + 
            gamma * (abs(U[1]) * interp3 - abs(U[0]) * interp4));
}

float convection_vT(float V[2], float T[5]) {
    float result = 0.0;
    
    float interp1 = (T[1] + T[3]) / 2;
    float interp2 = (T[4] + T[1]) / 2;
    float interp3 = (T[1] - T[3]) / 2;
    float interp4 = (T[4] - T[1]) /2;

    return inv_dy * ((V[1] * interp1 - V[0] * interp2) + 
            gamma * (abs(V[1]) * interp3 - abs(V[0]) * interp4));
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

float laplacian_5(float ar[5]) {
    float inv_dx2 = inv_dx * inv_dx;
    float inv_dy2 = inv_dy * inv_dy;
    float result = (ar[0] - 2. * ar[1] + ar[2]) * inv_dx2 +
                   (ar[3] - 2. * ar[1] + ar[4]) * inv_dy2 ;
    return result;
}