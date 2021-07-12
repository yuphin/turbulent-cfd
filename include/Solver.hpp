#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Boundary.hpp"
#include "Discretization.hpp"
#include "Domain.hpp"
#include "Fields.hpp"
#include "Grid.hpp"
#include "Utilities.hpp"


/**
 * @brief Parent class for all solver options (CPU, Cuda, Vulkan)
 *
 */
struct Solver {
    /**
     * @brief Calculation before solving the pressure equation
     *
     * Calculate dt, temperature, fluxes and right hand side for ppe in this order
     * Also apply necessary boundary conditions
     * 
     * @param[in] dt timestep size
     */
    virtual void solve_pre_pressure(Real &dt) = 0;
    virtual void solve_pressure(Real &res, uint32_t &it) = 0;
    virtual void solve_post_pressure() = 0;
    virtual void initialize() {}
    Solver() = default;

    /// Field object to be used
    Fields _field;
    /// Grid object to work on
    Grid _grid;
    /// Discretization object with necessary schemes
    Discretization _discretization;
    /// Vector of all boundaries of the problem
    std::vector<std::unique_ptr<Boundary>> _boundaries;

    /// Solver convergence tolerance
    Real _tolerance;

    /// 
    Real _omega;

    /// Precondition option for pcg solver
    int _preconditioner;

    /// Maximum number of iterations for the solver
    uint32_t _max_iter;

    // Turbulence modeling method
    int _turb_model = 0;
    
    Real _EPSIN;
    Real _KIN;
    
    friend class Logger;
    /// logger to output simulation information at runtime
    Logger logger = Logger();
    /// mpi parameters struct
    Params params;
    /// Solver type
    SolverType solver_type;
};
