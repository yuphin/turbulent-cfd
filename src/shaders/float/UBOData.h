
layout(push_constant) uniform UBOData {
	int imax;
    int jmax;
    int size;
    float nu;
    float alpha;
    float beta;
    float gx;
    float gy;
    float dx;
    float dy;
    float dx2;
    float dy2;
    float inv_dx;
    float inv_dy;
    float gamma;
    float PI;
    float UI;
    float VI;
    float tau;
    int num_diags;
    int num_fluid_cells;
    float wk;
    float weps;
};
