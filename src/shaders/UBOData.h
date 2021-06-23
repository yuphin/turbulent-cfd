layout(push_constant) uniform UBOData {
	int imax;
    int jmax;
    int size;
    float nu;
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
    int num_wgs;
    int num_diags;
};
