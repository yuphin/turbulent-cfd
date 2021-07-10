
layout(push_constant) uniform UBOData {
	int imax;
    int jmax;
    int size;
    double nu;
    double alpha;
    double beta;
    double gx;
    double gy;
    double dx;
    double dy;
    double dx2;
    double dy2;
    double inv_dx;
    double inv_dy;
    double gamma;
    double PI;
    double UI;
    double VI;
    double tau;
    int num_diags;
    int num_fluid_cells;
    double wk;
    double weps;
};

double dexp(double x) {
    double s = 1.0;
    for(int i = 10; i > 0; i--){
        s = 1 + x * s / i;
    }
    return s;
}

double dtanh(double x){
    double exp2x = dexp( 2 *x);
    return (exp2x - 1) / (exp2x + 1);
}
