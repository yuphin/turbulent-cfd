#ifndef __CUDACC__
#define __CUDACC__
#endif
#include "CUDASolver.hpp"
#include "Discretization.cu"
#include "device_launch_parameters.h"
#include <cuda.h>
#include <cuda_runtime.h>

#define BLK_SIZE 128
#define BLK_SIZE_2D 32

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
    if (code != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}
#define chk(ans)                                                                                                       \
    { gpuAssert((ans), __FILE__, __LINE__); }

int get_num_blks(int size) { return (size + BLK_SIZE - 1) / BLK_SIZE; }
dim3 get_num_blks_2d(int size_x, int size_y) {
    return dim3((size_x + BLK_SIZE_2D - 1) / BLK_SIZE_2D, (size_y + BLK_SIZE_2D - 1) / BLK_SIZE_2D);
}

template <typename T> inline void malloc_assign(T *dev_ptr, T val) {
    chk(cudaMalloc(&dev_ptr, sizeof(T)));
    chk(cudaMemcpy(dev_ptr, &val, 1, cudaMemcpyHostToDevice));
}

__global__ void init(Real *a, Real val, int size) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) {
        a[i] = val;
    }
}

__global__ void saxpy(Real *a, Real *x, Real *y, int size) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) {

        y[i] = *a * x[i] + y[i];
    }
}

__global__ void smaxpy(Real *a, Real *x, Real *y, int size) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) {
        y[i] = -(*a) * x[i] + y[i];
    }
}

__global__ void saxpy2(Real *a, Real *x, Real *y, int size) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) {
        y[i] = x[i] + *a * y[i];
    }
}

__global__ void vec_dot_vec(Real *a, Real *b, Real *o, int size) {
    __shared__ Real sdata[BLK_SIZE];
    uint32_t tid = threadIdx.x;
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i == 0) {
        *o = 0;
    }
    if (i < size) {
        sdata[tid] = a[i] * b[i];
    } else {
        sdata[tid] = 0;
    }
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] = sdata[tid] + sdata[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0) {
        atomicAdd(o, sdata[0]);
    }
}
// https://www.nvidia.com/docs/IO/66889/nvr-2008-004.pdf
__global__ void spmv_dia(Real *data, int *offsets, int num_rows, int num_cols, int num_diags, Real *x, Real *y) {
    int row = blockDim.x * blockIdx.x + threadIdx.x;
    if (row < num_rows) {
        Real dot = 0;
        y[row] = 0;
        for (int n = 0; n < num_diags; n++) {
            int col = row + offsets[n];
            Real val = data[num_rows * n + row];
            if (col >= 0 && col < num_cols) dot += val * x[col];
        }
        y[row] += dot;
    }
}

// See https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
__global__ void reduce_abs_max(Real *input, Real *output, int size) {
    extern __shared__ Real sdata[];
    uint32_t tid = threadIdx.x;
    uint32_t i = blockIdx.x * blockDim.x * 2 + threadIdx.x;
    Real curr_max = (i < size) ? real_abs(input[i]) : 0;
    if (i + blockDim.x < size) {
        curr_max = real_max(curr_max, real_abs(input[i + blockDim.x]));
    }
    sdata[tid] = curr_max;
    __syncthreads();
    for (uint32_t s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] = real_max(sdata[tid], sdata[tid + s]);
        }
        __syncthreads();
    }

    if (tid == 0) {
        output[blockIdx.x] = sdata[0];
    }
}

__global__ void reduce_min(Real *input, Real *output, int size) {
    extern __shared__ Real sdata[];
    uint32_t tid = threadIdx.x;
    uint32_t i = blockIdx.x * blockDim.x * 2 + threadIdx.x;
    Real curr_max = (i < size) ? input[i] : 0;
    if (i + blockDim.x < size) {
        curr_max = real_min(curr_max, (input[i + blockDim.x]));
    }
    sdata[tid] = curr_max;
    __syncthreads();
    for (uint32_t s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] = real_min(sdata[tid], sdata[tid + s]);
        }
        __syncthreads();
    }

    if (tid == 0) {
        output[blockIdx.x] = sdata[0];
    }
}

__global__ void reduce_max(Real *input, Real *output, int size) {
    extern __shared__ Real sdata[];
    uint32_t tid = threadIdx.x;
    uint32_t i = blockIdx.x * blockDim.x * 2 + threadIdx.x;
    Real curr_max = (i < size) ? input[i] : 0;
    if (i + blockDim.x < size) {
        curr_max = real_max(curr_max, (input[i + blockDim.x]));
    }
    sdata[tid] = curr_max;
    __syncthreads();
    for (uint32_t s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] = real_max(sdata[tid], sdata[tid + s]);
        }
        __syncthreads();
    }

    if (tid == 0) {
        output[blockIdx.x] = sdata[0];
    }
}

__global__ void scalar_div(Real *num, Real *denom, Real *o) { *o = *num / *denom; }
__global__ void scalar_cpy(Real *dst, Real *src) { *dst = *src; }

void solve_pcg(Real *A, int *A_offsets, int num_diag, Real *x, Real *b, Real *q, Real *d, Real *r, Real *r_dot_r_old,
               Real *r_dot_r, Real *z, Real *cg_beta, Real &delta_new, Real *cg_alpha, Real *d_dot_q, int precondition,
               Real *M, int *M_offsets, int m_num_diag, uint32_t &it, uint32_t max_iter, Real eps, int vec_size) {
    int num_blocks = get_num_blks(vec_size);
    cudaMemcpy(r, b, vec_size * sizeof(Real), cudaMemcpyDeviceToDevice);
    cudaMemset(x, 0, vec_size * sizeof(Real));
    if (precondition != -1) {
        spmv_dia<<<num_blocks, BLK_SIZE>>>(M, M_offsets, vec_size, vec_size, m_num_diag, r, d);
        vec_dot_vec<<<num_blocks, BLK_SIZE>>>(r, d, r_dot_r, vec_size);
    } else {
        cudaMemcpy(d, b, vec_size * sizeof(Real), cudaMemcpyDeviceToDevice);
        vec_dot_vec<<<num_blocks, BLK_SIZE>>>(r, r, r_dot_r, vec_size);
    }
    cudaMemcpy(&delta_new, r_dot_r, sizeof(Real), cudaMemcpyDeviceToHost);
    Real cond = delta_new * eps * eps;
    it = 0;
    while (it < max_iter && delta_new > cond) {
        // q <- A * d
        spmv_dia<<<num_blocks, BLK_SIZE>>>(A, A_offsets, vec_size, vec_size, num_diag, d, q);
        vec_dot_vec<<<num_blocks, BLK_SIZE>>>(d, q, d_dot_q, vec_size);
        // cg_alpha <- r_dot_r / d_dot_q
        scalar_div<<<1, 1>>>(r_dot_r, d_dot_q, cg_alpha);
        // x <- x + cg_alpha * d
        saxpy<<<num_blocks, BLK_SIZE>>>(cg_alpha, d, x, vec_size);
        // r <- r - cg_alpha * q
        smaxpy<<<num_blocks, BLK_SIZE>>>(cg_alpha, q, r, vec_size);
        scalar_cpy<<<1, 1>>>(r_dot_r_old, r_dot_r);
        if (precondition != -1) {
            // z <- M * r
            spmv_dia<<<num_blocks, BLK_SIZE>>>(M, M_offsets, vec_size, vec_size, m_num_diag, r, z);
            vec_dot_vec<<<num_blocks, BLK_SIZE>>>(r, z, r_dot_r, vec_size);
        } else {
            vec_dot_vec<<<num_blocks, BLK_SIZE>>>(r, r, r_dot_r, vec_size);
        }
        // cg_beta <- r_dot_r / r_dot_r_old
        scalar_div<<<1, 1>>>(r_dot_r, r_dot_r_old, cg_beta);
        if (precondition != -1) {
            // d <- z + cg_beta * d
            saxpy2<<<num_blocks, BLK_SIZE>>>(cg_beta, z, d, vec_size);
        } else {
            // d <- r + cg_beta *d
            saxpy2<<<num_blocks, BLK_SIZE>>>(cg_beta, r, d, vec_size);
        }
        it++;
        cudaMemcpy(&delta_new, r_dot_r, sizeof(Real), cudaMemcpyDeviceToHost);
    }
}


void get_abs_max(Real *input, Real *res, Real &out, int size) {
    int num_blks(get_num_blks(size));
    int smemsize = min(BLK_SIZE, size);
    reduce_abs_max<<<num_blks, BLK_SIZE, smemsize * sizeof(Real)>>>(input, res, size);
    while (num_blks != 1) {
        size = ceil(size / Real(BLK_SIZE));
        smemsize = min(BLK_SIZE, size);
        num_blks = get_num_blks(size);
        reduce_abs_max<<<num_blks, BLK_SIZE, smemsize * sizeof(Real)>>>(res, res, size);
    }
    cudaMemcpy(&out, res, sizeof(Real), cudaMemcpyDeviceToHost);
}

void solve_pcg2(Real *A, int *A_offsets, int num_diag, Real *x, Real *b, Real *q, Real *d, Real *r, Real *rho_old,
               Real *rho, Real *z, Real *cg_beta, Real *residual, Real &delta_new, Real *cg_alpha, Real *d_dot_q, int precondition,
               Real *M, int *M_offsets, int m_num_diag, uint32_t &it, uint32_t max_iter, Real eps, int vec_size) {
    int num_blocks = get_num_blks(vec_size);
    Real residual_out;
    cudaMemcpy(r, b, vec_size * sizeof(Real), cudaMemcpyDeviceToDevice);
    cudaMemset(x, 0, vec_size * sizeof(Real));
    get_abs_max(r, residual, residual_out, vec_size);
    if (residual_out == 0) {
        it = 0;
        return;
    }
    if (precondition != -1) {
        spmv_dia<<<num_blocks, BLK_SIZE>>>(M, M_offsets, vec_size, vec_size, m_num_diag, r, d);
        vec_dot_vec<<<num_blocks, BLK_SIZE>>>(r, d, rho, vec_size);
    } else {
        cudaMemcpy(d, b, vec_size * sizeof(Real), cudaMemcpyDeviceToDevice);
        vec_dot_vec<<<num_blocks, BLK_SIZE>>>(r, r, rho, vec_size);
    }
    Real residual_0 = residual_out;
    it = 0;
    while (it < max_iter) {
        // q <- A * d
        cudaMemset(q, 0, vec_size * sizeof(Real));
        spmv_dia<<<num_blocks, BLK_SIZE>>>(A, A_offsets, vec_size, vec_size, num_diag, d, q);
        vec_dot_vec<<<num_blocks, BLK_SIZE>>>(d, q, d_dot_q, vec_size);
        // cg_alpha <- rho / d_dot_q
        scalar_div<<<1, 1>>>(rho, d_dot_q, cg_alpha);
        // x <- x - cg_alpha * d
        smaxpy<<<num_blocks, BLK_SIZE>>>(cg_alpha, d, x, vec_size);
        // r <- r - cg_alpha * q
        smaxpy<<<num_blocks, BLK_SIZE>>>(cg_alpha, q, r, vec_size);
        get_abs_max(r, residual, residual_out, vec_size);
        if (residual_out <= eps) {
            break;
        }
        scalar_cpy<<<1, 1>>>(rho_old, rho);
        if (precondition != -1) {
            // z <- M * r
            spmv_dia<<<num_blocks, BLK_SIZE>>>(M, M_offsets, vec_size, vec_size, m_num_diag, r, z);
            vec_dot_vec<<<num_blocks, BLK_SIZE>>>(r, z, rho, vec_size);
        } else {
            vec_dot_vec<<<num_blocks, BLK_SIZE>>>(r, r, rho, vec_size);
        }
        // cg_beta <- rho / rho_old
        scalar_div<<<1, 1>>>(rho, rho_old, cg_beta);
        if (precondition != -1) {
            // d <- z + cg_beta * d
            saxpy2<<<num_blocks, BLK_SIZE>>>(cg_beta, z, d, vec_size);
        } else {
            // d <- r + cg_beta *d
            saxpy2<<<num_blocks, BLK_SIZE>>>(cg_beta, r, d, vec_size);
        }
        it++;
    }
    delta_new = residual_out / residual_0;
}

__global__ void sor_iter(Real *P, Real *RS, Real coeff, int *cell_type, int imax, int jmax, Real omega, Real inv_dx,
                         Real inv_dy, int parity) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    int is_fluid = at(cell_type, i, j);
    if (i >= imax || j >= jmax || is_fluid == 0) {
        return;
    }

    if (parity == 0 && ((i + j) % 2) == 0) {
        return;
    } else if (parity == 1 && ((i + j) % 2) == 1) {
        return;
    }
    Real p_stencil[4] = {at(P, i + 1, j), at(P, i - 1, j), at(P, i, j + 1), at(P, i, j - 1)};
    at(P, i, j) = (1 - omega) * at(P, i, j) + coeff * (sor_helper(p_stencil, inv_dx, inv_dy) - at(RS, i, j));
}

__global__ void calc_residual(Real *P, Real *RS, int *cell_type, int imax, int jmax, Real inv_dx, Real inv_dy,
                              Real *residual) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    int is_fluid = at(cell_type, i, j);
    if (i >= imax || j >= jmax || is_fluid == 0) {
        return;
    }

    Real rloc = 0;
    Real p_laplacian[5] = {at(P, i + 1, j), at(P, i, j), at(P, i - 1, j), at(P, i, j + 1), at(P, i, j - 1)};

    Real val = laplacian_5(p_laplacian, inv_dx, inv_dy) - at(RS, i, j);
    rloc += val * val;
    at(residual, i, j) = rloc;
}

__global__ void reduce_residual(Real *residual, Real *o, int size) {
    __shared__ Real sdata[BLK_SIZE];
    uint32_t tid = threadIdx.x;
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i == 0) {
        *o = 0;
    }
    if (i < size) {
        sdata[tid] = residual[i];
    } else {
        sdata[tid] = 0;
    }
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] = sdata[tid] + sdata[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0) {
        atomicAdd(o, sdata[0]);
    }
}

__global__ void negate_p(Real *p, int size) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) {
        p[i] = -p[i];
    }
}

__device__ Real vel_kernel(Real fg, Real dt, Real p[2], Real inv_dxy) { return fg - dt * inv_dxy * (p[1] - p[0]); }

__global__ void calc_vel(Real *u, Real *v, Real *p, Real *f, Real *g, int *cell_type, Real dt, int imax, int jmax,
                         Real dx, Real dy) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    int is_fluid = at(cell_type, i, j);
    if (i >= imax || j >= jmax || is_fluid == 0) {
        return;
    }
    Real inv_dx = 1 / dx;
    Real inv_dy = 1 / dy;

    Real p_diff_u[2] = {at(p, i, j), at(p, i + 1, j)};
    Real p_diff_v[2] = {at(p, i, j), at(p, i, j + 1)};
    /*  at(u, i, j) = vel_kernel(at(f, i, j), *dt, p_diff_u, inv_dx);
      at(v, i, j) = vel_kernel(at(g, i, j), *dt, p_diff_v, inv_dy);*/
    at(u, i, j) = at(f, i, j) - dt * inv_dx * (p_diff_u[1] - p_diff_u[0]);
    at(v, i, j) = at(g, i, j) - dt * inv_dy * (p_diff_v[1] - p_diff_v[0]);
}

__global__ void enforce_boundary(Real *u, int *row_start, int *col_idx, Real *mat, Real *rhs_vec, int size) {
    uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row < size) {
        Real sum = 0;
        for (int j = row_start[row]; j < row_start[row + 1]; j++) {
            sum += mat[j] * u[col_idx[j]];
        }
        u[row] = sum + 2 * rhs_vec[row];
    }
}

void uv_boundary(Real *u, Real *v, int *row_start_u, int *row_start_v, int *col_idx_u, int *col_idx_v, Real *mat_u,
                 Real *mat_v, Real *rhs_vec_u, Real *rhs_vec_v, int size) {
    int num_blks(get_num_blks(size));
    enforce_boundary<<<num_blks, BLK_SIZE>>>(u, row_start_u, col_idx_u, mat_u, rhs_vec_u, size);
    enforce_boundary<<<num_blks, BLK_SIZE>>>(v, row_start_v, col_idx_v, mat_v, rhs_vec_v, size);
}



void t_boundary(Real *t, int *row_start_t, int *col_idx_t, Real *mat_t, Real *rhs_vec_t, int size) {
    int num_blks(get_num_blks(size));
    enforce_boundary<<<num_blks, BLK_SIZE>>>(t, row_start_t, col_idx_t, mat_t, rhs_vec_t, size);
}

__global__ void calc_t(Real *u, Real *v, Real dx, Real dy, Real *t_new, Real *t_old, int *cell_type, Real alpha,
                       Real dt, Real gamma, int imax, int jmax) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    int is_fluid = at(cell_type, i, j);
    if (i >= imax || j >= jmax || is_fluid == 0) {
        return;
    }
    Real inv_dx = 1 / dx;
    Real inv_dy = 1 / dy;
    Real inv_dx2 = inv_dx * inv_dx;
    Real inv_dy2 = inv_dy * inv_dy;
    Real u_stencil[2] = {at(u, i - 1, j), at(u, i, j)};
    Real v_stencil[2] = {at(v, i, j - 1), at(v, i, j)};

    Real t_laplacian[5] = {at(t_old, i + 1, j), at(t_old, i, j), at(t_old, i - 1, j), at(t_old, i, j + 1),
                           at(t_old, i, j - 1)};

    at(t_new, i, j) = at(t_new, i, j) + dt * (alpha * laplacian_5(t_laplacian, inv_dx, inv_dy) -
                                              convection_uT(u_stencil, t_laplacian, inv_dx, gamma) -
                                              convection_vT(v_stencil, t_laplacian, inv_dy, gamma));
}

__global__ void calculate_nu_t(Real *NU_T, Real *K, Real *EPS, int *cell_type, Real _nu, int imax, int jmax) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    int is_fluid = at(cell_type, i, j);
    if (i >= imax || j >= jmax || is_fluid == 0) {
        return;
    }
    Real kij = at(K, i, j);
    Real epsij = at(EPS, i, j);
    at(NU_T, i, j) = 0.09 * kij * kij / epsij + _nu;
}

__global__ void calculate_nu_ij(Real *NU_I, Real *NU_J, Real *K, Real *EPS, int *cell_type, Real _nu, int imax,
                                int jmax) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    int is_fluid = at(cell_type, i, j);
    if (i >= imax || j >= jmax || is_fluid == 0) {
        return;
    }
    Real num_i = (at(K, i, j) + at(K, i + 1, j)) / 2;
    Real denom_i = (at(EPS, i, j) + at(EPS, i + 1, j)) / 2;

    Real num_j = (at(K, i, j) + at(K, i, j + 1)) / 2;
    Real denom_j = (at(EPS, i, j) + at(EPS, i, j + 1)) / 2;

    at(NU_I, i, j) = 0.09 * num_i * num_i / denom_i;
    at(NU_J, i, j) = 0.09 * num_j * num_j / denom_j;
}

__global__ void calculate_k_and_epsilon(Real *K_old, Real *EPS_old, Real *K, Real *EPS, Real *NU_T, Real *NU_I,
                                        Real *NU_J, Real *U, Real *V, int *cell_type, Real _nu, int imax, int jmax,
                                        Real dt, Real inv_dx, Real inv_dy) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    int is_fluid = at(cell_type, i, j);
    if (i >= imax || j >= jmax || is_fluid == 0) {
        return;
    }
    Real f2_coeff = 1;
    auto nut = at(NU_T, i, j);
    auto kij = at(K_old, i, j);
    auto eij = at(EPS_old, i, j);
    Real K_stencil[5] = {at(K_old, i + 1, j), at(K_old, i, j), at(K_old, i - 1, j), at(K_old, i, j + 1),
                         at(K_old, i, j - 1)};
    Real EPS_stencil[5] = {at(EPS_old, i + 1, j), at(EPS_old, i, j), at(EPS_old, i - 1, j), at(EPS_old, i, j + 1),
                           at(EPS_old, i, j - 1)};
    Real U_diff[2] = {at(U, i - 1, j), at(U, i, j)};
    Real V_diff[2] = {at(V, i, j - 1), at(V, i, j)};
    Real NU_I_diff[2] = {at(NU_I, i, j), at(NU_I, i - 1, j)};
    Real NU_J_diff[2] = {at(NU_J, i, j), at(NU_J, i, j - 1)};
    Real U_stencil[6] = {at(U, i, j),         at(U, i - 1, j), at(U, i, j + 1),
                         at(U, i - 1, j + 1), at(U, i, j - 1), at(U, i - 1, j - 1)};
    Real V_stencil[6] = {at(V, i, j),         at(V, i, j - 1), at(V, i + 1, j),
                         at(V, i + 1, j - 1), at(V, i - 1, j), at(V, i - 1, j - 1)};
    auto k1_1 = convecton_uKEPS(U_diff, K_stencil, inv_dx);
    auto k1_2 = convecton_vKEPS(V_diff, K_stencil, inv_dy);
    auto e1_1 = convecton_uKEPS(U_diff, EPS_stencil, inv_dx);
    auto e1_2 = convecton_vKEPS(V_diff, EPS_stencil, inv_dy);

    auto k2 = laplacian_nu(K_stencil, NU_I_diff, NU_J_diff, inv_dx, inv_dy, _nu, 1);
    auto e2 = laplacian_nu(EPS_stencil, NU_I_diff, NU_J_diff, inv_dx, inv_dy, _nu, 1.3);

    auto k3 = nut * mean_strain_rate_squared(U_stencil, V_stencil, inv_dx, inv_dy);
    auto e3 = 1.44 * eij * k3 / kij;
    auto e4 = f2_coeff * 1.92 * eij * eij / kij;
    auto kij_new = kij + dt * (-(k1_1 + k1_2) + k2 + k3 - eij);
    auto epsij_new = eij + dt * (-(e1_1 + e1_2) + e2 + e3 - e4);
    at(K, i, j) = kij_new;
    at(EPS, i, j) = epsij_new;
}

__global__ void calc_fg(Real *f, Real *g, Real *u, Real *v, bool calc_temp, Real dx, Real dy, Real *t, int *cell_type,
                        Real dt, Real gamma, Real nu, Real beta, Real gx, Real gy, int imax, int jmax) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    int is_fluid = at(cell_type, i, j);
    if (i >= imax || j >= jmax || is_fluid == 0) {
        return;
    }
    Real inv_dx = 1 / dx;
    Real inv_dy = 1 / dy;
    Real inv_dx2 = inv_dx * inv_dx;
    Real inv_dy2 = inv_dy * inv_dy;
    // 5-point + 1 stencil for U and V
    Real u_stencil[6] = {at(u, i + 1, j), at(u, i, j),     at(u, i - 1, j),
                         at(u, i, j + 1), at(u, i, j - 1), at(u, i - 1, j + 1)};
    Real v_stencil[6] = {at(v, i + 1, j), at(v, i, j),     at(v, i - 1, j),
                         at(v, i, j + 1), at(v, i, j - 1), at(v, i + 1, j - 1)};

    // Calculate fluxes
    at(f, i, j) = at(u, i, j) + dt * (nu * laplacian(u_stencil, inv_dx, inv_dy) -
                                      convection_u(u_stencil, v_stencil, inv_dx, inv_dy, gamma));
    at(g, i, j) = at(v, i, j) + dt * (nu * laplacian(v_stencil, inv_dx, inv_dy) -
                                      convection_v(u_stencil, v_stencil, inv_dx, inv_dy, gamma));

    if (calc_temp) {
        Real term1 = at(t, i, j) + at(t, i + 1, j);
        Real term2 = at(t, i, j) + at(t, i, j + 1);
        at(f, i, j) -= beta * dt / 2 * (term1)*gx;
        at(g, i, j) -= beta * dt / 2 * (term2)*gy;
    }
}

__global__ void calc_fg_turbulent(Real *f, Real *g, Real *u, Real *v, Real *NU_T, bool calc_temp, Real dx, Real dy,
                                  Real *t, int *cell_type, Real dt, Real gamma, Real nu, Real beta, Real gx, Real gy,
                                  int imax, int jmax) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    int is_fluid = at(cell_type, i, j);
    if (i >= imax || j >= jmax || is_fluid == 0) {
        return;
    }
    Real inv_dx = 1 / dx;
    Real inv_dy = 1 / dy;
    Real inv_dx2 = inv_dx * inv_dx;
    Real inv_dy2 = inv_dy * inv_dy;
    // 5-point + 1 stencil for U and V
    Real u_stencil[6] = {at(u, i + 1, j), at(u, i, j),     at(u, i - 1, j),
                         at(u, i, j + 1), at(u, i, j - 1), at(u, i - 1, j + 1)};
    Real v_stencil[6] = {at(v, i + 1, j), at(v, i, j),     at(v, i - 1, j),
                         at(v, i, j + 1), at(v, i, j - 1), at(v, i + 1, j - 1)};

    Real nu_term1 = (at(NU_T, i, j) + at(NU_T, i + 1, j)) / 2;
    Real nu_term2 = (at(NU_T, i, j) + at(NU_T, i, j + 1)) / 2;
    // Calculate fluxes
    at(f, i, j) = at(u, i, j) + dt * (nu_term1 * laplacian(u_stencil, inv_dx, inv_dy) -
                                      convection_u(u_stencil, v_stencil, inv_dx, inv_dy, gamma));
    at(g, i, j) = at(v, i, j) + dt * (nu_term2 * laplacian(v_stencil, inv_dx, inv_dy) -
                                      convection_v(u_stencil, v_stencil, inv_dx, inv_dy, gamma));

    if (calc_temp) {
        Real term1 = at(t, i, j) + at(t, i + 1, j);
        Real term2 = at(t, i, j) + at(t, i, j + 1);
        at(f, i, j) -= beta * dt / 2 * (term1)*gx;
        at(g, i, j) -= beta * dt / 2 * (term2)*gy;
    }
}

__global__ void fg_boundary(Real *f, Real *g, Real *u, Real *v, int imax, int jmax, uint32_t *neighborhood,
                            int *cell_type) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    int is_fluid = at(cell_type, i, j);
    if (i >= imax || j >= jmax || is_fluid == 1) {
        return;
    }
    uint32_t type = at(neighborhood, i, j) >> 8;
    uint32_t neighbors = at(neighborhood, i, j) & 0xFF;

    if ((neighbors & 0x1) == 1) {
        at(f, i, j) = at(u, i, j);
    }
    if ((neighbors & 0x2) == 2) {
        at(f, i - 1, j) = at(u, i - 1, j);
    }
    if ((neighbors & 0x4) == 4) {
        at(g, i, j) = at(v, i, j);
    }
    if ((neighbors & 0x8) == 8) {
        at(g, i, j - 1) = at(v, i, j - 1);
    }
}

__global__ void p_boundary(Real *p, int imax, int jmax, uint32_t *neighborhood, int *cell_type, Real PI) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    int is_fluid = at(cell_type, i, j);
    if (i >= imax || j >= jmax || is_fluid == 1) {
        return;
    }
    uint32_t type = at(neighborhood, i, j) >> 8;
    uint32_t neighbors = at(neighborhood, i, j) & 0xFF;
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
        if (!diag) {
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

__global__ void nu_t_boundary(Real *NU_T, Real *K, Real *EPS, int imax, int jmax, uint32_t *neighborhood,
                              int *cell_type, Real wk, Real weps, Real _nu) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    int is_fluid = at(cell_type, i, j);
    if (i >= imax || j >= jmax || is_fluid == 1) {
        return;
    }
    uint32_t type = at(neighborhood, i, j) >> 8;
    uint32_t neighbors = at(neighborhood, i, j) & 0xFF;
    Real k = 0;
    Real eps = 0;
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
        if (valid) {
            at(K, i, j) = k;
            at(EPS, i, j) = eps;
            at(NU_T, i, j) = 0.09 * k * k / eps + _nu;
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
        if (!diag) {
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
        if (diag || valid) {
            at(K, i, j) = k;
            at(EPS, i, j) = eps;
            at(NU_T, i, j) = 0.09 * k * k / eps + _nu;
        }
    }
}

void solve_sor(Real *P, Real *P_tmp, Real *P_residual, Real *P_residual_out, uint32_t *neighborhood, int imax, int jmax,
               Real *RS, int *cell_type, uint32_t &it, uint32_t max_iter, Real dx, Real dy, Real PI, Real tolerance,
               Real &res, int num_fluid_cells) {
    it = 0;
    const Real omega = 1.7;
    auto grid_size = imax * jmax;
    Real coeff = omega / (2 * (1 / (dx * dx) + 1 / (dy * dy)));
    dim3 blk_size_2d(BLK_SIZE_2D, BLK_SIZE_2D);
    dim3 num_blks_2d = get_num_blks_2d(imax, jmax);
    int num_blks_1d(get_num_blks(grid_size));
    Real inv_dx = 1 / dx;
    Real inv_dy = 1 / dy;
    res = REAL_MAX;
    while (it < max_iter && res > tolerance) {
        sor_iter<<<num_blks_2d, blk_size_2d>>>(P, RS, coeff, cell_type, imax, jmax, omega, inv_dx, inv_dy, 0);
        sor_iter<<<num_blks_2d, blk_size_2d>>>(P, RS, coeff, cell_type, imax, jmax, omega, inv_dx, inv_dy, 1);
        cudaMemset(P_residual, 0, grid_size * sizeof(Real));
        calc_residual<<<num_blks_2d, blk_size_2d>>>(P, RS, cell_type, imax, jmax, inv_dx, inv_dy, P_residual);
        reduce_residual<<<num_blks_1d, BLK_SIZE>>>(P_residual, P_residual_out, grid_size);
        p_boundary<<<num_blks_2d, blk_size_2d>>>(P, imax, jmax, neighborhood, cell_type, PI);
        cudaMemcpy(&res, P_residual_out, sizeof(Real), cudaMemcpyDeviceToHost);
        res = std::sqrt(res / num_fluid_cells);
        it++;
    }
}

__global__ void calc_rs(Real *f, Real *g, Real *rs, Real dx, Real dy, int imax, int jmax, Real dt, int *cell_type) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t j = blockIdx.y * blockDim.y + threadIdx.y;
    int is_fluid = at(cell_type, i, j);

    if (i >= imax || j >= jmax || is_fluid == 0) {
        return;
    }

    Real inv_dx = 1 / dx;
    Real inv_dy = 1 / dy;
    Real f_diff[2] = {at(f, i, j), at(f, i - 1, j)};
    Real g_diff[2] = {at(g, i, j), at(g, i, j - 1)};

    Real df = inv_dx * (f_diff[0] - f_diff[1]);
    Real dg = inv_dy * (g_diff[0] - g_diff[1]);
    at(rs, i, j) = (df + dg) * 1 / dt;
}

Real calculate_dt(int imax, int jmax, Real *u, Real *v, Real *u_residual, Real *v_residual, Real *nu_residual,
                  Real *k_residual, Real *eps_residual, Real *nu_t, Real *k, Real *eps, Real dx, Real dy, Real tau,
                  Real nu, Real alpha, bool calc_temp, bool turbulent) {
    // Calculate uv max
    int size = imax * jmax;
    int num_blks(get_num_blks(size));
    Real u_max_abs = 0;
    Real v_max_abs = 0;
    Real nu_min = REAL_MAX;
    Real k_max;
    Real eps_max;

    Real dx2 = dx * dx;
    Real dy2 = dy * dy;
    int smemsize = min(BLK_SIZE, size);
    std::vector<Real> ucpu(imax * jmax);
    std::vector<Real> vcpu(imax * jmax);
    std::vector<Real> ures(imax * jmax);
    std::vector<Real> vres(imax * jmax);

    reduce_abs_max<<<num_blks, BLK_SIZE, smemsize * sizeof(Real)>>>(u, u_residual, size);
    reduce_abs_max<<<num_blks, BLK_SIZE, smemsize * sizeof(Real)>>>(v, v_residual, size);
    if (turbulent) {
        reduce_min<<<num_blks, BLK_SIZE, smemsize * sizeof(Real)>>>(nu_t, nu_residual, size);
        reduce_max<<<num_blks, BLK_SIZE, smemsize * sizeof(Real)>>>(k, k_residual, size);
        reduce_max<<<num_blks, BLK_SIZE, smemsize * sizeof(Real)>>>(eps, eps_residual, size);
    }
    while (num_blks != 1) {
        size = ceil(size / Real(BLK_SIZE));
        smemsize = min(BLK_SIZE, size);
        num_blks = get_num_blks(size);
        reduce_abs_max<<<num_blks, BLK_SIZE, smemsize * sizeof(Real)>>>(u_residual, u_residual, size);
        reduce_abs_max<<<num_blks, BLK_SIZE, smemsize * sizeof(Real)>>>(v_residual, v_residual, size);
        if (turbulent) {
            reduce_min<<<num_blks, BLK_SIZE, smemsize * sizeof(Real)>>>(nu_t, nu_residual, size);
            reduce_max<<<num_blks, BLK_SIZE, smemsize * sizeof(Real)>>>(k, k_residual, size);
            reduce_max<<<num_blks, BLK_SIZE, smemsize * sizeof(Real)>>>(eps, eps_residual, size);
        }
    }
    cudaMemcpy(&u_max_abs, u_residual, sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(&v_max_abs, v_residual, sizeof(Real), cudaMemcpyDeviceToHost);
    if (turbulent) {
        cudaMemcpy(&nu_min, nu_residual, sizeof(Real), cudaMemcpyDeviceToHost);
        cudaMemcpy(&k_max, k_residual, sizeof(Real), cudaMemcpyDeviceToHost);
        cudaMemcpy(&eps_max, eps_residual, sizeof(Real), cudaMemcpyDeviceToHost);
    }
    Real min_cond = std::min(dx / u_max_abs, dy / v_max_abs);

    nu_min = (nu_min == REAL_MAX || nu_min == 0) ? nu : nu_min;

    if (nu_min != 0) {
        Real cond_spatial = 1.0 / (2.0 * nu) * ((dx2 * dy2) / (dx2 + dy2));
        min_cond = std::min(min_cond, cond_spatial);
    }

    if (calc_temp) {
        Real inv_dx = 1 / dx;
        Real inv_dx2 = inv_dx * inv_dx;
        Real inv_dy = 1 / dy;
        Real inv_dy2 = inv_dy * inv_dy;
        Real cond_temp = 1 / (2 * alpha * (inv_dx2 + inv_dy2));
        min_cond = std::min(min_cond, cond_temp);
    }
    if (turbulent) {
        Real cond_5 = 1 / (2 * k_max * (1 / dx2 + 1 / dy2));
        Real cond_6 = 1 / (2 * eps_max * (1 / dx2 + 1 / dy2));
        min_cond = std::min(min_cond, cond_5);
        min_cond = std::min(min_cond, cond_6);
    }

    return tau * min_cond;
}
void CudaSolver::initialize() {
    auto grid_x = _grid.imaxb();
    auto grid_y = _grid.jmaxb();
    auto grid_size = grid_x * grid_y;
    build_pcg_matrix(_field, _grid, _boundaries, A_pcg, U_pcg, V_pcg, T_pcg, U_RHS, V_RHS, T_RHS, U_fixed, V_fixed,
                     T_fixed);
    // Preprocess
    std::vector<int> is_fluid(_grid.imaxb() * _grid.jmaxb(), 0);
    for (const auto &current_cell : _grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();
        is_fluid[_grid.imaxb() * j + i] = 1;
    }

    std::vector<BoundaryData> neighbors(_grid.imaxb() * _grid.jmaxb());
    for (const auto &boundary : _boundaries) {
        uint32_t type = boundary->get_type();
        auto cells = boundary->_cells;
        for (auto &cell : *cells) {
            int i = cell->i();
            int j = cell->j();
            BoundaryData data;
            uint32_t type = boundary->get_type();
            data.neighborhood |= type << 8;
            // data.idx = _grid.imaxb() * j + i;
            if (cell->is_border(border_position::RIGHT)) {
                data.neighborhood |= 1;
            }
            if (cell->is_border(border_position::LEFT)) {
                data.neighborhood |= 2;
            }
            if (cell->is_border(border_position::TOP)) {
                data.neighborhood |= 4;
            }
            if (cell->is_border(border_position::BOTTOM)) {
                data.neighborhood |= 8;
            }
            if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::TOP)) {
                data.neighborhood |= 16;
            }
            if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::BOTTOM)) {
                data.neighborhood |= 32;
            }
            if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::TOP)) {
                data.neighborhood |= 64;
            }
            if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::BOTTOM)) {
                data.neighborhood |= 128;
            }
            neighbors[j * _grid.imaxb() + i] = data;
        }
    }
    DiagonalSparseMatrix<Real> A_matrix_diag =
        create_diagonal_matrix(A_pcg, _grid.imaxb(), _grid.jmaxb(), {-_grid.imaxb(), -1, 0, 1, _grid.imaxb()});
    DiagonalSparseMatrix<Real> A_precond_diag;
    if (_preconditioner != -1) {
        A_precond_diag = create_preconditioner_spai(A_pcg, _grid, _preconditioner);
    }
    num_offsets_a = A_matrix_diag.offsets.size();
    num_offsets_m = A_precond_diag.offsets.size();
    auto t_matrix_data = T_fixed.value.data();
    auto t_matrix_size = T_fixed.value.size();
    auto t_row_start_data = T_fixed.rowstart.data();
    auto t_row_start_size = T_fixed.rowstart.size();
    auto t_col_idx_data = T_fixed.colindex.data();
    auto t_col_idx_size = T_fixed.colindex.size();
    auto t_rhs_data = T_RHS.data();
    auto t_rhs_size = T_RHS.size();
    auto u_matrix_data = U_fixed.value.data();
    auto u_matrix_size = U_fixed.value.size();
    auto v_matrix_data = V_fixed.value.data();
    auto v_matrix_size = V_fixed.value.size();
    auto u_row_start_data = U_fixed.rowstart.data();
    auto u_row_start_size = U_fixed.rowstart.size();
    auto u_col_idx_data = U_fixed.colindex.data();
    auto u_col_idx_size = U_fixed.colindex.size();
    auto v_row_start_data = V_fixed.rowstart.data();
    auto v_row_start_size = V_fixed.rowstart.size();
    auto v_col_idx_data = V_fixed.colindex.data();
    auto v_col_idx_size = V_fixed.colindex.size();
    auto u_rhs_data = U_RHS.data();
    auto u_rhs_size = U_RHS.size();
    auto v_rhs_data = V_RHS.data();
    auto v_rhs_size = V_RHS.size();

    cudaMalloc(&U, grid_size * sizeof(Real));
    cudaMalloc(&V, grid_size * sizeof(Real));
    cudaMalloc(&F, grid_size * sizeof(Real));
    cudaMalloc(&G, grid_size * sizeof(Real));
    cudaMalloc(&P, grid_size * sizeof(Real));
    cudaMalloc(&P_temp, grid_size * sizeof(Real));

    cudaMalloc(&RS, grid_size * sizeof(Real));
    cudaMalloc(&U_residual, grid_size * sizeof(Real));
    cudaMalloc(&V_residual, grid_size * sizeof(Real));
    if (_turb_model != 0) {
        cudaMalloc(&NU_residual, grid_size * sizeof(Real));
        cudaMalloc(&K_residual, grid_size * sizeof(Real));
        cudaMalloc(&EPS_residual, grid_size * sizeof(Real));
        cudaMalloc(&NU_T, grid_size * sizeof(Real));
        cudaMalloc(&NU_I, grid_size * sizeof(Real));
        cudaMalloc(&NU_J, grid_size * sizeof(Real));
        cudaMalloc(&K, grid_size * sizeof(Real));
        cudaMalloc(&K_old, grid_size * sizeof(Real));
        cudaMalloc(&EPS, grid_size * sizeof(Real));
        cudaMalloc(&EPS_old, grid_size * sizeof(Real));

        cudaMemset(NU_residual, 0, grid_size * sizeof(Real));
        cudaMemset(K_residual, 0, grid_size * sizeof(Real));
        cudaMemset(EPS_residual, 0, grid_size * sizeof(Real));
        cudaMemset(NU_T, 0, grid_size * sizeof(Real));
        cudaMemset(NU_I, 0, grid_size * sizeof(Real));
        cudaMemset(NU_J, 0, grid_size * sizeof(Real));
        cudaMemcpy(K, _field._K._container.data(), grid_size * sizeof(Real), cudaMemcpyHostToDevice);
        cudaMemcpy(EPS, _field._EPS._container.data(), grid_size * sizeof(Real), cudaMemcpyHostToDevice);
    }
    cudaMalloc(&P_residual, grid_size * sizeof(Real));
    cudaMalloc(&cell_type, grid_size * sizeof(int));
    cudaMalloc(&row_start_u, u_row_start_size * sizeof(int));
    cudaMalloc(&row_start_v, v_row_start_size * sizeof(int));
    cudaMalloc(&row_start_t, t_row_start_size * sizeof(int));
    cudaMalloc(&col_idx_u, u_col_idx_size * sizeof(int));
    cudaMalloc(&col_idx_v, v_col_idx_size * sizeof(int));
    cudaMalloc(&mat_u, u_matrix_size * sizeof(Real));
    cudaMalloc(&mat_v, v_matrix_size * sizeof(Real));
    cudaMalloc(&rhs_vec_u, u_rhs_size * sizeof(Real));
    cudaMalloc(&rhs_vec_v, v_rhs_size * sizeof(Real));
    cudaMalloc(&neighborhood, neighbors.size() * sizeof(uint32_t));

    cudaMalloc(&A, A_matrix_diag.data.size() * sizeof(Real));
    cudaMalloc(&A_offsets, A_matrix_diag.offsets.size() * sizeof(uint32_t));
    if (_preconditioner != -1) {
        cudaMalloc(&M, A_precond_diag.data.size() * sizeof(Real));
        cudaMalloc(&M_offsets, A_precond_diag.offsets.size() * sizeof(uint32_t));
    }
    cudaMalloc(&q, grid_size * sizeof(Real));
    cudaMalloc(&d, grid_size * sizeof(Real));
    cudaMalloc(&r, grid_size * sizeof(Real));
    cudaMalloc(&z, grid_size * sizeof(Real));
    cudaMalloc(&r_dot_r, sizeof(Real));
    cudaMalloc(&r_dot_r_old, sizeof(Real));
    cudaMalloc(&d_dot_q, sizeof(Real));
    cudaMalloc(&p_residual_out, sizeof(Real));
    cudaMalloc(&cg_alpha, sizeof(Real));
    cudaMalloc(&cg_beta, sizeof(Real));

    chk(cudaMemcpy(U, _field._U._container.data(), grid_size * sizeof(Real), cudaMemcpyHostToDevice));
    chk(cudaMemcpy(V, _field._V._container.data(), grid_size * sizeof(Real), cudaMemcpyHostToDevice));
    chk(cudaMemcpy(P, _field._P._container.data(), grid_size * sizeof(Real), cudaMemcpyHostToDevice));
    chk(cudaMemcpy(mat_u, u_matrix_data, u_matrix_size * sizeof(Real), cudaMemcpyHostToDevice));
    chk(cudaMemcpy(row_start_u, u_row_start_data, u_row_start_size * sizeof(int), cudaMemcpyHostToDevice));
    chk(cudaMemcpy(col_idx_u, u_col_idx_data, u_col_idx_size * sizeof(int), cudaMemcpyHostToDevice));
    chk(cudaMemcpy(rhs_vec_u, u_rhs_data, u_rhs_size * sizeof(Real), cudaMemcpyHostToDevice));
    chk(cudaMemcpy(mat_v, v_matrix_data, v_matrix_size * sizeof(Real), cudaMemcpyHostToDevice));
    chk(cudaMemcpy(row_start_v, v_row_start_data, v_row_start_size * sizeof(int), cudaMemcpyHostToDevice));
    chk(cudaMemcpy(col_idx_v, v_col_idx_data, v_col_idx_size * sizeof(int), cudaMemcpyHostToDevice));
    chk(cudaMemcpy(rhs_vec_v, v_rhs_data, v_rhs_size * sizeof(Real), cudaMemcpyHostToDevice));

    chk(cudaMemcpy(neighborhood, neighbors.data(), neighbors.size() * sizeof(uint32_t), cudaMemcpyHostToDevice));

    chk(cudaMemcpy(A, A_matrix_diag.data.data(), A_matrix_diag.data.size() * sizeof(Real), cudaMemcpyHostToDevice));
    chk(cudaMemcpy(A_offsets, A_matrix_diag.offsets.data(), A_matrix_diag.offsets.size() * sizeof(int),
                   cudaMemcpyHostToDevice));
    if (_field.calc_temp) {
        cudaMalloc(&T, grid_size * sizeof(Real));
        cudaMalloc(&T_temp, grid_size * sizeof(Real));
        cudaMalloc(&mat_t, t_matrix_size * sizeof(Real));
        cudaMalloc(&rhs_vec_t, t_rhs_size * sizeof(Real));
        cudaMalloc(&col_idx_t, t_col_idx_size * sizeof(int));

        cudaMemcpy(mat_t, t_matrix_data, t_matrix_size * sizeof(Real), cudaMemcpyHostToDevice);
        cudaMemcpy(row_start_t, t_row_start_data, t_row_start_size * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(col_idx_t, t_col_idx_data, t_col_idx_size * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(rhs_vec_t, t_rhs_data, t_rhs_size * sizeof(Real), cudaMemcpyHostToDevice);
        chk(cudaMemcpy(T, _field._T._container.data(), grid_size * sizeof(Real), cudaMemcpyHostToDevice));
    }

    if (_preconditioner != -1) {
        chk(cudaMemcpy(M, A_precond_diag.data.data(), A_precond_diag.data.size() * sizeof(Real),
                       cudaMemcpyHostToDevice));
        chk(cudaMemcpy(M_offsets, A_precond_diag.offsets.data(), A_precond_diag.offsets.size() * sizeof(int),
                       cudaMemcpyHostToDevice));
    }
    chk(cudaMemcpy(cell_type, is_fluid.data(), is_fluid.size() * sizeof(int), cudaMemcpyHostToDevice));
}
void CudaSolver::solve_pre_pressure(Real &dt) {
    auto grid_x = _grid.imaxb();
    auto grid_y = _grid.jmaxb();
    auto grid_size = grid_x * grid_y;
    dim3 num_blks_1d(get_num_blks(grid_size));
    dim3 blk_size_2d(BLK_SIZE_2D, BLK_SIZE_2D);
    dim3 num_blks_2d = get_num_blks_2d(grid_x, grid_y);
    dt = calculate_dt(_grid.imaxb(), _grid.jmaxb(), U, V, U_residual, V_residual, NU_residual, K_residual, EPS_residual,
                      NU_T, K, EPS, _grid.dx(), _grid.dy(), _field._tau, _field._nu, _field._alpha, _field.calc_temp,
                      _turb_model != 0);
    _field._dt = dt;
    uv_boundary(U, V, row_start_u, row_start_v, col_idx_u, col_idx_v, mat_u, mat_v, rhs_vec_u, rhs_vec_v, grid_size);
    if (_field.calc_temp) {
        t_boundary(T, row_start_t, col_idx_t, mat_t, rhs_vec_t, grid_size);
        chk(cudaMemcpy(T_temp, T, grid_size * sizeof(Real), cudaMemcpyDeviceToDevice));
        calc_t<<<num_blks_2d, blk_size_2d>>>(U, V, _grid.dx(), _grid.dy(), T, T_temp, cell_type, _field._alpha, dt,
                                             _discretization._gamma, _grid.imaxb(), _grid.jmaxb());
    }
    std::vector<Real> fcpu(grid_size);
    std::vector<Real> gcpu(grid_size);
    if (_turb_model != 0) {
        calc_fg_turbulent<<<num_blks_2d, blk_size_2d>>>(F, G, U, V, NU_T, _field.calc_temp, _grid.dx(), _grid.dy(), T,
                                                        cell_type, dt, _discretization._gamma, _field._nu, _field._beta,
                                                        _field._gx, _field._gy, grid_x, grid_y);
    } else {
        calc_fg<<<num_blks_2d, blk_size_2d>>>(F, G, U, V, _field.calc_temp, _grid.dx(), _grid.dy(), T, cell_type, dt,
                                              _discretization._gamma, _field._nu, _field._beta, _field._gx, _field._gy,
                                              grid_x, grid_y);
    }
    fg_boundary<<<num_blks_2d, blk_size_2d>>>(F, G, U, V, grid_x, grid_y, neighborhood, cell_type);
    calc_rs<<<num_blks_2d, blk_size_2d>>>(F, G, RS, _grid.dx(), _grid.dy(), grid_x, grid_y, dt, cell_type);
}

void CudaSolver::solve_pressure(Real &res, uint32_t &it) {
    if (solver_type == SolverType::PCG) {
        auto grid_x = _grid.imaxb();
        auto grid_y = _grid.jmaxb();
        auto grid_size = grid_x * grid_y;
        int num_blks(get_num_blks(grid_size));
        solve_pcg2(A, A_offsets, num_offsets_a, P, RS, q, d, r, r_dot_r_old, r_dot_r, z, cg_beta, U_residual, res, cg_alpha, d_dot_q,
                  _preconditioner, M, M_offsets, num_offsets_m, it, _max_iter, _tolerance,
                  _grid.imaxb() * _grid.jmaxb());
    } else if (solver_type == SolverType::SOR) {
        solve_sor(P, P_temp, P_residual, p_residual_out, neighborhood, _grid.imaxb(), _grid.jmaxb(), RS, cell_type, it,
                  _max_iter, _grid.dx(), _grid.dy(), _field._PI, _tolerance, res, _grid.fluid_cells().size());
    }
}

void CudaSolver::solve_post_pressure() {
    auto grid_x = _grid.imaxb();
    auto grid_y = _grid.jmaxb();
    auto grid_size = grid_x * grid_y;
    int num_blks(get_num_blks(grid_size));
    dim3 blk_size_2d(BLK_SIZE_2D, BLK_SIZE_2D);
    dim3 num_blks_2d = get_num_blks_2d(grid_x, grid_y);
    calc_vel<<<num_blks_2d, blk_size_2d>>>(U, V, P, F, G, cell_type, _field._dt, grid_x, grid_y, _grid.dx(),
                                           _grid.dy());
    if (_turb_model == 1) {
        chk(cudaMemcpy(K_old, K, grid_size * sizeof(Real), cudaMemcpyDeviceToDevice));
        chk(cudaMemcpy(EPS_old, EPS, grid_size * sizeof(Real), cudaMemcpyDeviceToDevice));
        calculate_nu_t<<<num_blks_2d, blk_size_2d>>>(NU_T, K, EPS, cell_type, _field._nu, grid_x, grid_y);
        calculate_nu_ij<<<num_blks_2d, blk_size_2d>>>(NU_I, NU_J, K, EPS, cell_type, _field._nu, grid_x, grid_y);
        calculate_k_and_epsilon<<<num_blks_2d, blk_size_2d>>>(K_old, EPS_old, K, EPS, NU_T, NU_I, NU_J, U, V, cell_type,
                                                              _field._nu, grid_x, grid_y, _field._dt, 1 / _grid.dx(),
                                                              1 / _grid.dy());
        // TODO : Implement KIN and EPSIN
        nu_t_boundary<<<num_blks_2d, blk_size_2d>>>(NU_T, K, EPS, grid_x, grid_y, neighborhood, cell_type, 0, 0,
                                                    _field._nu);
    }
    chk(cudaMemcpy(_field._U._container.data(), U, grid_size * sizeof(Real), cudaMemcpyDeviceToHost));
    chk(cudaMemcpy(_field._V._container.data(), V, grid_size * sizeof(Real), cudaMemcpyDeviceToHost));
    chk(cudaMemcpy(_field._P._container.data(), P, grid_size * sizeof(Real), cudaMemcpyDeviceToHost));
    if (_field.calc_temp) {
        chk(cudaMemcpy(_field._T._container.data(), T, grid_size * sizeof(Real), cudaMemcpyDeviceToHost));
    }
}

CudaSolver::~CudaSolver() {
    cudaFree(U);
    cudaFree(V);
    cudaFree(F);
    cudaFree(G);
    cudaFree(P);
    cudaFree(T);
    cudaFree(T_temp);
    cudaFree(RS);
    cudaFree(U_residual);
    cudaFree(V_residual);
    if (_turb_model != 0) {
        cudaFree(NU_residual);
        cudaFree(K_residual);
        cudaFree(EPS_residual);
        cudaFree(NU_T);
        cudaFree(NU_I);
        cudaFree(NU_J);
        cudaFree(K);
        cudaFree(K_old);
        cudaFree(EPS_old);
        cudaFree(EPS);
    }
    cudaFree(P_residual);
    cudaFree(cell_type);
    cudaFree(P);
    cudaFree(P_temp);
    cudaFree(row_start_u);
    cudaFree(row_start_v);
    cudaFree(row_start_t);
    cudaFree(col_idx_u);
    cudaFree(col_idx_v);
    cudaFree(col_idx_t);
    cudaFree(mat_u);
    cudaFree(mat_v);
    cudaFree(mat_t);
    cudaFree(rhs_vec_u);
    cudaFree(rhs_vec_v);
    cudaFree(rhs_vec_t);
    cudaFree(neighborhood);
    cudaFree(A);
    cudaFree(A_offsets);
    if (_preconditioner != -1) {
        cudaFree(M);
        cudaFree(M_offsets);
    }
    cudaFree(q);
    cudaFree(d);
    cudaFree(r);
    cudaFree(z);
    cudaFree(r_dot_r);
    cudaFree(r_dot_r_old);
    cudaFree(d_dot_q);
    cudaFree(p_residual_out);
    cudaFree(cg_alpha);
    cudaFree(cg_beta);
}
