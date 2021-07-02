#include "CPUSolver.hpp"
#include "Communication.hpp"

void CPUSolver::solve_pre_pressure(Real &dt) {

    // Select dt
    dt = _field.calculate_dt(_grid, _calc_temp);

    // Enforce velocity boundary conditions
    for (auto &boundary : _boundaries) {
        boundary->enforce_uv(_field);
    }

    if (_calc_temp) {
        // Enforce temperature boundary conditions
        for (const auto &boundary : _boundaries) {
            boundary->enforce_t(_field);
        }
        // Compute temperatures
        _field.calculate_temperatures(_grid);

        // Communicate temperatures
        Communication::communicate(&params, _field.t_matrix());
    }

    // Compute F & G and enforce boundary conditions
    _field.calculate_fluxes(_grid, _calc_temp);
    for (const auto &boundary : _boundaries) {
        boundary->enforce_fg(_field);
    }
    // Communicate F and G
    Communication::communicate(&params, _field.f_matrix());
    Communication::communicate(&params, _field.g_matrix());
    // Set RHS of PPE
    _field.calculate_rs(_grid);
}

void CPUSolver::solve_pressure(Real &res, uint32_t &it) {
   
    // Perform pressure solve
    res = _pressure_solver->solve(_field, _grid, _boundaries, params, _max_iter, _tolerance, it);
}

void CPUSolver::solve_post_pressure() {
    // Compute u^(n+1) & v^(n+1)
    _field.calculate_velocities(_grid);
    // Communicate velocities
    Communication::communicate(&params, _field.u_matrix());
    Communication::communicate(&params, _field.v_matrix());
}
