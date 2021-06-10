#include "PressureSolver.hpp"

#include <cmath>
#include <iostream>

SOR::SOR(Real omega) : _omega(omega) {}

Real SOR::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries, uint32_t iters,
                Real tolerance) {

    Real dx = grid.dx();
    Real dy = grid.dy();

    Real coeff = _omega / (2.0f * (1.0 / (dx * dx) + 1.0 / (dy * dy))); // = _omega * h^2 / 4.0, if dx == dy == h

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        field.p(i, j) = (1.0 - _omega) * field.p(i, j) +
                        coeff * (Discretization::sor_helper(field.p_matrix(), i, j) - field.rs(i, j));
    }

    Real res = 0.0;
    Real rloc = 0.0;

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        Real val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
        rloc += (val * val);
    }
    {
        res = rloc / (grid.fluid_cells().size());
        res = std::sqrt(res);
    }

    return res;
}

PCG::PCG(int dim_x, int dim_y, Real dx, Real dy, Fields &field, Grid &grid,
         const std::vector<std::unique_ptr<Boundary>> &boundaries)
    : dim(dim_x * dim_y), dim_x(dim_x), dim_y(dim_y), A(dim) {
    build_matrix(A, dx, dy, field, grid, boundaries);
    create_diagonal_matrix();
}

Real PCG::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries, uint32_t iters,
                Real tolerance) {
    Real pcg_residual = 1e10;
    int pcg_iters = -1;
    static SparsePCGSolver<Real> solver;
    solver.set_solver_parameters(tolerance, iters);
    int dim_x = this->dim_x;
    auto at = [dim_x](int i, int j) { return j * dim_x + i; };
    Real inv_dx2 = 1 / (grid.dx() * grid.dx());
    Real inv_dy2 = 1 / (grid.dy() * grid.dy());
    Real div = inv_dx2 + inv_dy2;
    std::vector<Real> p(dim, 0);

    solver.solve(A, field._RS._container, p, pcg_residual, pcg_iters, 0);
    if (pcg_iters == 0) {
        int a = 4;
    }
    for (int i = 0; i < field._P.size(); i++) {
        field._P._container[i] = -p[i];
    }
    std::cout << pcg_iters;
    // field._P._container = p;
    return pcg_residual;
}

void PCG::build_matrix(SparseMatrix<Real> &A, Real dx, Real dy, Fields &field, Grid &grid,
                       const std::vector<std::unique_ptr<Boundary>> &boundaries) {
    int dim_x = this->dim_x;
    Real inv_dx2 = 1 / (dx * dx);
    Real inv_dy2 = 1 / (dy * dy);
    Real div = inv_dx2 + inv_dy2;
    auto at = [dim_x](int i, int j) { return j * dim_x + i; };
    /* A.set_element(0, 0, 1);
     A.set_element(at(dim_x - 1, dim_y - 1), at(dim_x - 1, dim_y - 1), 1);
     for (int j = 0; j < dim_y; j++) {
         for (int i = 0; i < dim_x; i++) {
             auto loc = at(i, j);
             A.set_element(loc, loc, 1);
         }
     }*/
    for (auto &boundary : boundaries) {
        for (auto &cell : *boundary->_cells) {
            int i = cell->i();
            int j = cell->j();
            auto loc = at(i, j);
            if (boundary->get_type() == 0) {
                // Outlet
                A.set_element(loc, loc, 1);
                field.rs(i, j) = field._PI;

            } else if (boundary->get_type() == 1 || boundary->get_type() == 2) {
                field.rs(i, j) = 0;
                A.set_element(loc, loc, 1 * inv_dx2);
                // Inlet & noslip
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
                    A.set_element(loc, loc, 0.5 * (dim));
                    A.set_element(loc, at(i + 1, j), -0.5 * inv_dx2);
                    A.set_element(loc, at(i, j + 1), -0.5 * inv_dy2);
                }
                if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::BOTTOM)) {
                    A.set_element(loc, loc, 0.5 * (dim));
                    A.set_element(loc, at(i + 1, j), -0.5 * inv_dx2);
                    A.set_element(loc, at(i, j - 1), -0.5 * inv_dy2);
                }
                if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::TOP)) {
                    A.set_element(loc, loc, 0.5 * (dim));
                    A.set_element(loc, at(i - 1, j), -0.5 * inv_dx2);
                    A.set_element(loc, at(i, j + 1), -0.5 * inv_dy2);
                }
                if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::BOTTOM)) {
                    A.set_element(loc, loc, 0.5 * (dim));
                    A.set_element(loc, at(i - 1, j), -0.5 * inv_dx2);
                    A.set_element(loc, at(i, j - 1), -0.5 * inv_dy2);
                }
            }
        }
    }

    for (auto current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();
        auto loc = at(i, j);
        A.set_element(loc, loc, 2 * div);
        A.set_element(loc, at(i - 1, j), -inv_dx2);
        A.set_element(loc, at(i + 1, j), -inv_dx2);
        A.set_element(loc, at(i, j + 1), -inv_dy2);
        A.set_element(loc, at(i, j - 1), -inv_dy2);
    }
}

void PCG::create_diagonal_matrix() {
    A_diag.dim = dim;
    A_diag.offsets = {-dim_x, -1, 0, 1, dim_x};
    A_diag.data.resize(5 * dim);
    int cnt = 0;
    for (auto r = 0; r < dim; r++) {
        for (int n = 0; n < A_diag.num_diags; n++) {
            auto c = A_diag.offsets[n] + r;
            if (c >= 0 && c < dim) {
                A_diag.data[dim * n + r] = A(r, c);
            }
        }
    }
}
