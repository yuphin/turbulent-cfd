#include "CPUSolver.hpp"
#include "Communication.hpp"

void CPUSolver::initialize() {
    if (solver_type == SolverType::PCG) {
        build_pcg_matrix();
    }
}

void CPUSolver::solve_pre_pressure(Real &dt) {

    // Select dt
    dt = _field.calculate_dt(_grid, _field.calc_temp, _turb_model);

    // Enforce velocity boundary conditions
    for (auto &boundary : _boundaries) {
        boundary->enforce_uv(_field);
    }

    if (_field.calc_temp) {
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
    _field.calculate_fluxes(_grid, _field.calc_temp, _turb_model != 0);
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
    res = REAL_MAX;
    it = 0;
    if (solver_type == SolverType::SOR) {
        while (it < _max_iter && res > _tolerance) {
            res = solve_sor();
            // Enforce boundary conditions
            for (const auto &boundary : _boundaries) {
                boundary->enforce_p(_field);
            }
            it++;
        }

    } else if (solver_type == SolverType::PCG) {
        res = solve_pcg(it); 
    }
}

void CPUSolver::solve_post_pressure() {
    // Compute u^(n+1) & v^(n+1)
    _field.calculate_velocities(_grid);
    // Communicate velocities
    Communication::communicate(&params, _field.u_matrix());
    Communication::communicate(&params, _field.v_matrix());

     if (_turb_model != 0) {
        // Compute turbulent viscosity and set boundary conditions
        _field.calculate_nu_t(_grid, _turb_model);
        for (const auto &boundary : _boundaries) {
            boundary->enforce_nu_t(_field, _turb_model, _grid.dx(), _grid.dy());
        }
        // Communicate turbulence quantities
        Communication::communicate(&params, _field.nu_t_matrix());
        Communication::communicate(&params, _field._NU_I);
        Communication::communicate(&params, _field._NU_J);
        Communication::communicate(&params, _field.k_matrix());
        Communication::communicate(&params, _field.eps_matrix());
    }
}

Real CPUSolver::solve_pcg(uint32_t &it) {
    Real pcg_residual = 1e10;
    int pcg_iters = -1;
    static SparsePCGSolver<Real> solver;
    solver.set_solver_parameters(_tolerance, _max_iter);
    int dim_x = _grid.imaxb();
    int dim = _grid.domain().total_size;
    auto at = [dim_x](int i, int j) { return j * dim_x + i; };
    Real inv_dx2 = 1 / (_grid.dx() * _grid.dx());
    Real inv_dy2 = 1 / (_grid.dy() * _grid.dy());
    Real div = inv_dx2 + inv_dy2;
    solver.solve(A, _field._RS._container, _field._P._container, pcg_residual, pcg_iters, 0);
    it = pcg_iters;
    return pcg_residual;

}

Real CPUSolver::solve_sor() {
    Real dx = _grid.dx();
    Real dy = _grid.dy();

    Real coeff = _omega / (2.0f * (1.0 / (dx * dx) + 1.0 / (dy * dy))); // = _omega * h^2 / 4.0, if dx == dy == h

    for (auto currentCell : _grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        _field.p(i, j) = (1.0 - _omega) * _field.p(i, j) +
                         coeff * (Discretization::sor_helper(_field.p_matrix(), i, j) - _field.rs(i, j));
    }
    Real res = 0.0;
    Real rloc = 0.0;

    for (auto currentCell : _grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        Real val = Discretization::laplacian(_field.p_matrix(), i, j) - _field.rs(i, j);
        rloc += (val * val);
    }
    // Exchange pressure
    Communication::communicate(&params, _field.p_matrix());
    // Compute global residual
    int global_cells;
    res = Communication::reduce_all(rloc, MPI_SUM);
    global_cells = Communication::reduce_all(_grid.fluid_cells().size(), MPI_SUM);
    res = res / global_cells;
    res = std::sqrt(res);
    return res;
}

void CPUSolver::build_pcg_matrix() {
    int dim_x = _grid.imaxb();
    auto dx = _grid.dx();
    auto dy = _grid.dy();
    Real inv_dx2 = 1 / (dx * dx);
    Real inv_dy2 = 1 / (dy * dy);
    Real div = inv_dx2 + inv_dy2;
    auto at = [dim_x](int i, int j) { return j * dim_x + i; };
    A.resize(_grid.domain().total_size);
    for (auto current_cell : _grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();
        auto loc = at(i, j);
        A.set_element(loc, loc, 2 * div);
        A.set_element(loc, at(i - 1, j), -inv_dx2);
        A.set_element(loc, at(i + 1, j), -inv_dx2);
        A.set_element(loc, at(i, j + 1), -inv_dy2);
        A.set_element(loc, at(i, j - 1), -inv_dy2);
    }
    for (auto &boundary : _boundaries) {
        for (auto &cell : *boundary->_cells) {
            int i = cell->i();
            int j = cell->j();
            int id = cell->id();
            auto loc = at(i, j);
            switch (boundary->get_type()) {
            case 0: {
                // Outlet
                A.set_element(loc, loc, 1);
                _field.rs(i, j) = _field._PI;
            } break;
            case 1: {
                // Inlet
                auto inlet_vel_u = static_cast<InletBoundary *>(boundary.get())->_inlet_U[id];
                auto inlet_vel_v = static_cast<InletBoundary *>(boundary.get())->_inlet_V[id];
                _field.rs(i, j) = 0;
                A.set_element(loc, loc, 1 * inv_dx2);
                // Apply noslip for now
                if (cell->is_border(border_position::RIGHT)) {
                    A.set_element(loc, loc, inv_dx2);
                    A.set_element(loc, at(i + 1, j), -inv_dx2);

                } else if (cell->is_border(border_position::LEFT)) {
                    A.set_element(loc, loc, inv_dx2);
                    A.set_element(loc, at(i - 1, j), -inv_dx2);

                } else if (cell->is_border(border_position::TOP)) {
                    A.set_element(loc, loc, inv_dy2);
                    A.set_element(loc, at(i, j + 1), -inv_dy2);

                } else if (cell->is_border(border_position::BOTTOM)) {
                    A.set_element(loc, loc, inv_dy2);
                    A.set_element(loc, at(i, j - 1), -inv_dy2);
                }
                if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::TOP)) {
                    A.set_element(loc, loc, 0.5 * (div));
                    A.set_element(loc, at(i + 1, j), -0.5 * inv_dx2);
                    A.set_element(loc, at(i, j + 1), -0.5 * inv_dy2);
                }
                if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::BOTTOM)) {
                    A.set_element(loc, loc, 0.5 * (div));
                    A.set_element(loc, at(i + 1, j), -0.5 * inv_dx2);
                    A.set_element(loc, at(i, j - 1), -0.5 * inv_dy2);
                }
                if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::TOP)) {
                    A.set_element(loc, loc, 0.5 * (div));
                    A.set_element(loc, at(i - 1, j), -0.5 * inv_dx2);
                    A.set_element(loc, at(i, j + 1), -0.5 * inv_dy2);
                }
                if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::BOTTOM)) {
                    A.set_element(loc, loc, 0.5 * (div));
                    A.set_element(loc, at(i - 1, j), -0.5 * inv_dx2);
                    A.set_element(loc, at(i, j - 1), -0.5 * inv_dy2);
                }
            } break;
            case 2: {
                // NoSlip
                _field.rs(i, j) = 0;
                A.set_element(loc, loc, 1 * inv_dx2);
                auto wall_vel_u = static_cast<NoSlipWallBoundary *>(boundary.get())->_wall_velocity[id];
                auto wall_vel_v = static_cast<NoSlipWallBoundary *>(boundary.get())->_wall_velocity[id];
                // Apply noslip for now
                if (cell->is_border(border_position::RIGHT)) {
                    A.set_element(loc, loc, inv_dx2);
                    A.set_element(loc, at(i + 1, j), -inv_dx2);
                } else if (cell->is_border(border_position::LEFT)) {
                    // Pressure
                    A.set_element(loc, loc, inv_dx2);
                    A.set_element(loc, at(i - 1, j), -inv_dx2);
                } else if (cell->is_border(border_position::TOP)) {
                    // Pressure
                    A.set_element(loc, loc, inv_dy2);
                    A.set_element(loc, at(i, j + 1), -inv_dy2);
                } else if (cell->is_border(border_position::BOTTOM)) {
                    // Pressure
                    A.set_element(loc, loc, inv_dy2);
                    A.set_element(loc, at(i, j - 1), -inv_dy2);
                }
                if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::TOP)) {
                    // Pressure
                    A.set_element(loc, loc, 0.5 * (div));
                    A.set_element(loc, at(i + 1, j), -0.5 * inv_dx2);
                    A.set_element(loc, at(i, j + 1), -0.5 * inv_dy2);
                }
                if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::BOTTOM)) {
                    // Pressure
                    A.set_element(loc, loc, 0.5 * (div));
                    A.set_element(loc, at(i + 1, j), -0.5 * inv_dx2);
                    A.set_element(loc, at(i, j - 1), -0.5 * inv_dy2);
                }
                if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::TOP)) {
                    // Pressure
                    A.set_element(loc, loc, 0.5 * (div));
                    A.set_element(loc, at(i - 1, j), -0.5 * inv_dx2);
                    A.set_element(loc, at(i, j + 1), -0.5 * inv_dy2);
                }
                if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::BOTTOM)) {
                    // Pressure
                    A.set_element(loc, loc, 0.5 * (div));
                    A.set_element(loc, at(i - 1, j), -0.5 * inv_dx2);
                    A.set_element(loc, at(i, j - 1), -0.5 * inv_dy2);
                }
            } break;
            case 3: {
                // FreeSlip
            } break;
            default:
                break;
            }
        }
    }
}
