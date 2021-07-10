#pragma once
#include "Solver.hpp"
#include "PCGUtils.h"

#if USE_FLOATS
#define real_max(v1, v2) fmaxf(v1, v2)
#define real_min(v1, v2) fminf(v1, v2)
#define real_abs(v1) fabsf(v1)
#define real_tanh(v1) tanhf(v1)
#define real_pow(v1, v2) powf(v1, v2)
#define real_sqrt(v1) sqrtf(v1)
#else
#define real_max(v1, v2) fmax(v1,v2)
#define real_min(v1, v2) fmin(v1,v2)
#define real_abs(v1) fabs(v1)
#define real_tanh(v1) tanh(v1)
#define real_pow(v1, v2) pow(v1, v2)
#define real_sqrt(v1) sqrt(v1)
#endif

struct CudaSolver : public Solver {
    CudaSolver() = default;
    ~CudaSolver();
    void initialize() override;
    virtual void solve_pre_pressure(Real &dt) override;
    virtual void solve_pressure(Real &res, uint32_t &it) override;
    virtual void solve_post_pressure() override;

  private:
    Real *U;
    Real *V;
    Real *F;
    Real *G;
    Real *P;
    Real *P_temp;
    Real *T;
    Real *T_temp;
    Real *RS;
    Real *U_residual;
    Real *V_residual;
    Real *NU_residual;
    Real *K_residual;
    Real *EPS_residual;
    Real *P_residual;
    int *cell_type;
    int *row_start_u;
    int *row_start_v;
    int *row_start_t;
    int *col_idx_u;
    int *col_idx_v;
    int *col_idx_t;
    Real *mat_u;
    Real *mat_v;
    Real *mat_t;
    Real *rhs_vec_u;
    Real *rhs_vec_v;
    Real *rhs_vec_t;
    uint32_t *neighborhood;

    Real *A;
    int *A_offsets;
    Real *M;
    int *M_offsets;
    Real *q;
    Real *d;
    Real *r;
    Real *z;
    Real *r_dot_r;
    Real *r_dot_r_old;
    Real *d_dot_q;
    Real *p_residual_out;
    Real *cg_alpha;
    Real *cg_beta;
    int num_offsets_a;
    int num_offsets_m;

    Real *NU_T;
    Real *NU_I;
    Real *NU_J;
    Real *K;
    Real *K_old;
    Real *EPS;
    Real *EPS_old;
    Real *dists;
    Real *S;

    Real *dx;
    Real *dy;
    int *imax;
    int *jmax;
    Real *tau;
    Real *nu;
    Real *alpha;
    Real *gamma;
    int *calc_temp;

    SparseMatrix<Real> A_pcg;
    SparseMatrix<Real> U_pcg;
    SparseMatrix<Real> V_pcg;
    SparseMatrix<Real> T_pcg;
    FixedSparseMatrix<Real> U_fixed;
    FixedSparseMatrix<Real> V_fixed;
    FixedSparseMatrix<Real> T_fixed;
    DiagonalSparseMatrix<Real> spai;
    std::vector<Real> U_RHS;
    std::vector<Real> V_RHS;
    std::vector<Real> T_RHS;
    

};