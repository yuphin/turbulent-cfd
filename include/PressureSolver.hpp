#pragma once

#include "Boundary.hpp"
#include <Utilities.hpp>
#include "Fields.hpp"
#include "Grid.hpp"
#include <lib/pcgsolver.h>
#include <utility>
/**
 * @brief Abstract class for pressure Poisson equation solver
 *
 */
class PressureSolver {
  public:
    PressureSolver() = default;
    virtual ~PressureSolver() = default;

    /**
     * @brief Solve the pressure equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid to be used
     * @param[in] boundary to be used
     */
    virtual Real solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries,
                       uint32_t iters, Real tolerance) = 0;
};

/**
 * @brief Successive Over-Relaxation algorithm for solution of pressure Poisson
 * equation
 *
 */
class SOR : public PressureSolver {
  public:
    SOR() = default;

    /**
     * @brief Constructor of SOR solver
     *
     * @param[in] relaxation factor
     */
    SOR(Real omega);

    virtual ~SOR() = default;

    /**
     * @brief Solve the pressure equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid to be used
     * @param[in] boundary to be used
     */
    virtual Real solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries,
                       uint32_t iters, Real tolerance);

  private:
    Real _omega;
};

class PCG : public PressureSolver {
  public:
    PCG(int dim_x, int dim_y, Real dx, Real dy, Fields &field, Grid &grid,
        const std::vector<std::unique_ptr<Boundary>> &boundaries);

    virtual Real solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries,
                       uint32_t iters, Real tolerance);
    DiagonalSparseMatrix<Real> A_diag;
    SparseMatrix<Real> &get_a() { return A; }
    FixedSparseMatrix<Real> U_fixed;
    FixedSparseMatrix<Real> V_fixed;
    std::vector<Real> U_RHS;
    std::vector<Real> V_RHS;
  private:
    int dim;
    int dim_x;
    int dim_y;
    SparseMatrix<Real> A;
    SparseMatrix<Real> U;
    SparseMatrix<Real> V;
    void build_matrix(Real dx, Real dy, Fields &field, Grid &grid,
                      const std::vector<std::unique_ptr<Boundary>> &boundaries);
    void create_diagonal_matrix();
};
