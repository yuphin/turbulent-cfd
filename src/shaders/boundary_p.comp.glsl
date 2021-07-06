#version 450
#extension GL_KHR_shader_subgroup_arithmetic : enable
#extension GL_GOOGLE_include_directive: require
#include "UBOData.h"
layout(local_size_x = 32, local_size_y = 32, local_size_z = 1) in;

#define at(AR, I, J) AR[imax * (J) + (I)]


layout(binding = 5) readonly buffer CellType
{
	float cell_type[];
};
layout(binding = 7) buffer P
{
	float p[];
};

layout(binding = 9) buffer Neighborhood {
    uint neighborhood[];
};


void main() {
	uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    float is_fluid = at(cell_type, i, j);
    if(i >= imax || j >= jmax || is_fluid == 1) {
        return;
    }
	uint type = at(neighborhood, i, j) >> 8;
	uint neighbors = at(neighborhood, i, j) & 0xFF;
	 if (type == 0) { // Outlet
        at(p, i, j) = PI;
    } else {
        int diag = 0;
        if ((neighbors & 0x10) == 16) { // Right + top
            diag = 1;
            at(p, i, j) = (at(p, i + 1, j) + at(p, i, j + 1)) / 2;
        }
        if ((neighbors & 0x20) == 32) { // Right + bottom
            diag = 1;
            at(p, i, j) = (at(p, i + 1, j) + at(p, i, j - 1)) / 2;
        }
        if ((neighbors & 0x40) == 64) { // Left + top
            diag = 1;
            at(p, i, j) = (at(p, i - 1, j) + at(p, i, j + 1)) / 2;
        }
        if ((neighbors & 0x80) == 128) { // Left + bottom
            diag = 1;
            at(p, i, j) = (at(p, i - 1, j) + at(p, i, j - 1)) / 2;
        }
        if (diag == 0) {
            if ((neighbors & 0x1) == 1) { // Right
                at(p, i, j) = at(p, i + 1, j);
            }
            if ((neighbors & 0x2) == 2) { // Left
                at(p, i, j) = at(p, i - 1, j);
            }
            if ((neighbors & 0x4) == 4) { // Top
                at(p, i, j) = at(p, i, j + 1);
            }
            if ((neighbors & 0x8) == 8) { // Bottom
                at(p, i, j) = at(p, i, j - 1);
            }
        }
    }
}