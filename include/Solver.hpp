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
struct Solver {
    virtual void solve_pre_pressure(Real &dt) = 0;
    virtual void solve_pressure(Real &res, uint32_t &it) = 0;
    virtual void solve_post_pressure() = 0;
    virtual void initialize() {}
    Solver() = default;
    Fields _field;
    Grid _grid;
    Discretization _discretization;
    std::vector<std::unique_ptr<Boundary>> _boundaries;

    /// Solver convergence tolerance
    Real _tolerance;

    Real _omega;

    int _preconditioner;

    /// Maximum number of iterations for the solver
    uint32_t _max_iter;

    // Turbulence modeling method
    int _turb_model = 0;

    
    friend class Logger;
    Logger logger = Logger();
    Params params;
    SolverType solver_type;
};
