#pragma once
#include "Solver.hpp"
#include <lib/pcgsolver.h>

/**
 * @brief Main logic to solve the problem on a CPU with MPI support
 *
 */
struct CPUSolver : public Solver {
    CPUSolver() = default;

    /**
     * @brief Do extra initialization if necessary
     */
    virtual void initialize() override;

    /**
     * @brief Calculation before solving the pressure equation
     *
     * Calculate dt, temperature, fluxes and right hand side for ppe in this order
     * Also apply necessary boundary conditions
     * 
     * @param[in] dt timestep size
     */
    virtual void solve_pre_pressure(Real &dt) override;

    /**
     * @brief Solve the pressure equation with specified numerical solver
     * 
     * @param[in] res Residual
     * @param[in] it Number of iterations of the solver
     */
    virtual void solve_pressure(Real &res, uint32_t  &it) override;

    /**
     * @brief Calculation after solving the pressure equation
     *
     * Calculate velocities and turbulence values
     */
    virtual void solve_post_pressure() override;

private:
    /**
     * @brief Solve the ppe using pcg algorithm
     * 
     * @param[in] it number of iterations
     * @param[out] res The global residual
     */
    Real solve_pcg(uint32_t &it);

    /**
     * @brief Solve the ppe using sor algorithm
     * 
     * @param[in] it number of iterations
     * @param[out] res The global residual
     */
    Real solve_sor();

    /// Build the pcg matrix
    void build_pcg_matrix();

    /// PCG matrix for the PCG solver
    SparseMatrix<Real> A;

};