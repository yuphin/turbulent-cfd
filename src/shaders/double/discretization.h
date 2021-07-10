#define at(AR, I, J) AR[imax * (J) + (I)]
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
double convection_u(double U[6], double V[6]) {
    double result = 0.0;
    double interp1 = (U[0] + U[1]) / 2;
    double interp2 = (U[2] + U[1]) / 2;
    double interp3 = interp1 * 2;
    double interp4 = (U[1] - U[0]) / 2;
    double interp5 = interp2 * 2;
    double interp6 = (U[2] - U[1]) /2;

    double interp7 = (V[0] + V[1]) / 2;
    double interp8 = (U[3] + U[1]) / 2;
    double interp9 = (V[5] + V[4]) / 2;
    double interp10 = (U[4] + U[1]) /2;
    double interp11 = interp7 * 2;
    double interp12 = (U[1] - U[3]) /2;
    double interp13 = interp9 * 2;
    double interp14 = (U[4] - U[1]) /2;
   
    // dU^2/dx
    double result_fd = (interp1 * interp1 - interp2 * interp2) * inv_dx;
    double result_dc = gamma * (abs(interp3) / 2 * interp4 - abs(interp5) /2 * interp6) * inv_dx; 
    result += result_fd + result_dc;

    // dUV/dy
    result_fd = (interp7 * interp8 - interp9 * interp10) * inv_dy;
    result_dc = gamma * (abs(interp11) / 2 * interp12 -
                        abs(interp13) / 2 * interp14 ) * inv_dy;
    result += result_fd + result_dc;

    return result;

}

double convection_uT(double U[2], double T[5]) {
    double result = 0.0;
    
    double interp1 = (T[1] + T[0]) / 2;
    double interp2 = (T[2] + T[1]) / 2;
    double interp3 = (T[1] - T[0]) / 2;
    double interp4 = (T[2] - T[1]) /2;

    return inv_dx * ((U[1] * interp1 - U[0] * interp2) + 
            gamma * (abs(U[1]) * interp3 - abs(U[0]) * interp4));
}

double convection_vT(double V[2], double T[5]) {
    double result = 0.0;
    
    double interp1 = (T[1] + T[3]) / 2;
    double interp2 = (T[4] + T[1]) / 2;
    double interp3 = (T[1] - T[3]) / 2;
    double interp4 = (T[4] - T[1]) /2;

    return inv_dy * ((V[1] * interp1 - V[0] * interp2) + 
            gamma * (abs(V[1]) * interp3 - abs(V[0]) * interp4));
}


double convection_v(double U[6], double V[6]) {
    double result = 0.0;
    double interp1 = (V[3] + V[1]) / 2;
    double interp2 = (V[4] + V[1]) / 2;
    double interp3 = interp1 * 2;
    double interp4 = (V[1] - V[3]) / 2;
    double interp5 = interp2 * 2;
    double interp6 = (V[4] - V[1]) /2;

    double interp7 = (U[3] + U[1]) / 2;
    double interp8 = (V[0] + V[1]) / 2;
    double interp9 = (U[5] + U[2]) / 2;
    double interp10 = (V[2] + V[1]) /2;
    double interp11 = interp7 * 2;
    double interp12 = (V[1] - V[0]) /2;
    double interp13 = interp9 * 2;
    double interp14 = (V[2] - V[1]) /2;
   
    // dU^2/dx
    double result_fd = (interp1 * interp1 - interp2 * interp2) * inv_dy;
    double result_dc = gamma * (abs(interp3) / 2 * interp4 - abs(interp5) * interp6) * inv_dy; 
    result += result_fd + result_dc;

    // dUV/dy
    result_fd = (interp7 * interp8 - interp9 * interp10) * inv_dx;
    result_dc = gamma * (abs(interp11) / 2 * interp12 -
                        abs(interp13) / 2 * interp14 ) * inv_dx;
    result += result_fd + result_dc;
    return result;
}

double laplacian(double ar[6]) {
    double inv_dx2 = inv_dx * inv_dx;
    double inv_dy2 = inv_dy * inv_dy;
    double result = (ar[0] - 2. * ar[1] + ar[2]) * inv_dx2 +
                   (ar[3] - 2. * ar[1] + ar[4]) * inv_dy2 ;
    return result;
}

double laplacian_5(double ar[5]) {
    double inv_dx2 = inv_dx * inv_dx;
    double inv_dy2 = inv_dy * inv_dy;
    double result = (ar[0] - 2. * ar[1] + ar[2]) * inv_dx2 +
                   (ar[3] - 2. * ar[1] + ar[4]) * inv_dy2 ;
    return result;
}