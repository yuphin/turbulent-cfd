#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;
layout(constant_id = 0) const uint TURB_MODEL = 0;
#define at(AR, I, J) AR[imax * (J) + (I)]


layout(binding = 5) readonly buffer CellType
{
	double cell_type[];
};

layout(binding = 9) buffer Neighborhood {
    uint neighborhood[];
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

void main() {
	uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    double is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 1) {
        return;
    }
	uint type = at(neighborhood, i, j) >> 8;
	uint neighbors = at(neighborhood, i, j) & 0xFF;
    double k = 0;
    double eps = 0;
	   if (type == 1) { // Inlet
        int valid = 0;
        if ((neighbors & 0x1) == 1) { // Right
            k = 2 * wk - at(K, i + 1, j);
            eps = 2 * weps - at(EPS, i + 1, j);
            valid = 1;
        }
        if ((neighbors & 0x2) == 2) { // Left
            k = 2 * wk - at(K, i - 1, j);
            eps = 2 * weps - at(EPS, i - 1, j);
            valid = 1;
        }
        if ((neighbors & 0x4) == 4) { // Top
            k = 2 * wk - at(K, i, j + 1);
            eps = 2 * weps - at(EPS, i, j + 1);
            valid = 1;
        }
        if ((neighbors & 0x8) == 8) { // Bottom
            k = 2 * wk - at(K, i, j - 1);
            eps = 2 * weps - at(EPS, i, j - 1);
            valid = 1;
        }
        if (valid == 1) {
            at(K, i, j) = k;
            at(EPS, i, j) = eps;
            // at(NU_T, i, j) = 0.09 * k * k / eps + _nu;
        }
    } else { // Other
        int diag = 0;
        int valid = 0;
        if ((neighbors & 0x10) == 16) { // Right + top
            diag = 1;
            k = (at(K, i + 1, j) + at(K, i, j + 1)) / 2;
            eps = (at(EPS, i + 1, j) + at(EPS, i, j + 1)) / 2;
        }
        if ((neighbors & 0x20) == 32) { // Right + bottom
            diag = 1;
            k = (at(K, i + 1, j) + at(K, i, j - 1)) / 2;
            eps = (at(EPS, i + 1, j) + at(EPS, i, j - 1)) / 2;
        }
        if ((neighbors & 0x40) == 64) { // Left + top
            diag = 1;
            k = (at(K, i - 1, j) + at(K, i, j + 1)) / 2;
            eps = (at(EPS, i - 1, j) + at(EPS, i, j + 1)) / 2;
        }
        if ((neighbors & 0x80) == 128) { // Left + bottom
            diag = 1;
            k = (at(K, i - 1, j) + at(K, i, j - 1)) / 2;
            eps = (at(EPS, i - 1, j) + at(EPS, i, j - 1)) / 2;
        }
        if (diag == 0) {
            if ((neighbors & 0x1) == 1) { // Right
                k = at(K, i + 1, j);
                eps = at(EPS, i + 1, j);
                valid = 1;
            }
            if ((neighbors & 0x2) == 2) { // Left
                k = at(K, i - 1, j);
                eps = at(EPS, i - 1, j);
                valid = 1;
            }
            if ((neighbors & 0x4) == 4) { // Top
                k = at(K, i, j + 1);
                eps = at(EPS, i, j + 1);
                valid = 1;
            }
            if ((neighbors & 0x8) == 8) { // Bottom
                k = at(K, i, j - 1);
                eps = at(EPS, i, j - 1);
                valid = 1;
            }
        }
        if (diag==1 || valid == 1) {
            at(K, i, j) = k;
            at(EPS, i, j) = eps;
            if (TURB_MODEL == 1) {
                at(NU_T, i, j) = 0.09 * k * k / eps + nu;
            } else if (TURB_MODEL == 2 || TURB_MODEL == 3) {
                at(NU_T, i, j) = k / eps + nu;
            }
        }
    }
}