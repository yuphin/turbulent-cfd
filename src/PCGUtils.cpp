#include "PCGUtils.h"
DiagonalSparseMatrix<Real> create_diagonal_matrix(const SparseMatrix<Real> &A, int dim_x, int dim_y,
                                                  const std::vector<int> offsets) {
    DiagonalSparseMatrix<Real> A_diag;
    A_diag.dim = A.n;
    A_diag.offsets = offsets;
    A_diag.num_diags = offsets.size();
    A_diag.data.resize(offsets.size() * A.n);
    int cnt = 0;
    for (auto r = 0; r < A.n; r++) {
        for (int n = 0; n < A_diag.num_diags; n++) {
            auto c = A_diag.offsets[n] + r;
            if (c >= 0 && c < A.n) {
                A_diag.data[A.n * n + r] = A(r, c);
            }
        }
    }
    return A_diag;
}

DiagonalSparseMatrix<Real> create_preconditioner_spai(const SparseMatrix<Real> &A, Grid &_grid, int precond_type) {
    std::vector<int> diag_offsets;
    auto dim_x = _grid.imaxb();
    auto dim_y = _grid.jmaxb();
    SparseMatrix<Real> result(_grid.domain().total_size);
    if (precond_type == 0) {
        const Real omg = 0.5;
        SparseMatrix<Real> diag(_grid.domain().total_size);
        SparseMatrix<Real> j_a(_grid.domain().total_size);
        SparseMatrix<Real> j_a_j(_grid.domain().total_size);
        // From "Algorithm for Sparse Approximate Inverse Preconditioners in the Conjugate Gradient Method",
        // https://interval.louisiana.edu/reliable-computing-journal/volume-19/reliable-computing-19-pp-120-126.pdf
        for (int i = 0; i < A.n; i++) {
            diag.set_element(i, i, 1 / A(i, i));
        }
        mat_mat_multiply(diag, A, j_a);
        mat_mat_multiply(j_a, diag, j_a_j);
        for (int i = 0; i < j_a_j.n; i++) {
            for (const int j : j_a_j.index[i]) {
                auto elem = j_a_j(i, j);
                Real diag = 0;
                if (i == j) {
                    diag = 2 / A(i, i);
                }
                result.set_element(i, j, diag - elem);
            }
        }
        for (int i = 0; i < result.n; i++) {
            for (const int j : result.index[i]) {
                auto a = result(i, j);
                Real d = 0;
                if (i == j) {
                    d = diag(i, i);
                }
                result.set_element(i, j, omg * d + (1 - omg) * result(i, j));
            }
        }
        diag_offsets = {-dim_x, -1, 0, 1, dim_x};
    } else if (precond_type == 1) {
        const Real omg = 0.5;
        SparseMatrix<Real> K_inv(dim_x * dim_y);
        SparseMatrix<Real> K_inv_T(dim_x * dim_y);
        // SSOR preconditioner from "Parallel preconditioned conjugate gradient algorithm on GPU"
        for (int i = 0; i < A.n; i++) {
            for (const int j : A.index[i]) {
                if (j > i) {
                    continue;
                }
                auto delta = i == j;
                auto term = omg * A(i, j) / A(j, j);
                K_inv.set_element(i, j, sqrt(2 - omg) * sqrt(omg / A(i, i)) * (delta - term));
            }
        }
        mat_transpose(K_inv, K_inv_T);
        // Calculate M_inv = K_inv_T * K_inv
        mat_mat_multiply(K_inv_T, K_inv, result);
        // diag_offsets = {-dim_x + 1, 0, dim_x - 1};
        diag_offsets = {-dim_x, -dim_x + 1, -1, 0, 1, dim_x - 1, dim_x};
    } else if (precond_type == 2) {
        for (int i = 0; i < A.n; i++) {
            result.set_element(i, i, 1 / A(i, i));
        }
        diag_offsets = {0};
    }
    return create_diagonal_matrix(result, dim_x, dim_y, diag_offsets);
}

void build_pcg_matrix(Fields &_field, Grid &_grid, const std::vector<std::unique_ptr<Boundary>> &_boundaries,
                      SparseMatrix<Real> &A, SparseMatrix<Real> &U, SparseMatrix<Real> &V, SparseMatrix<Real> &T,
                      std::vector<Real> &U_RHS, std::vector<Real> &V_RHS, std::vector<Real> &T_RHS,
                      FixedSparseMatrix<Real> &U_fixed, FixedSparseMatrix<Real> &V_fixed,
                      FixedSparseMatrix<Real> &T_fixed) {
    int dim_x = _grid.imaxb();
    auto dx = _grid.dx();
    auto dy = _grid.dy();
    Real inv_dx2 = 1 / (dx * dx);
    Real inv_dy2 = 1 / (dy * dy);
    Real div = inv_dx2 + inv_dy2;
    auto dim = _grid.domain().total_size;
    auto at = [dim_x](int i, int j) { return j * dim_x + i; };
    U_RHS.resize(dim);
    V_RHS.resize(dim);
    A.resize(dim);
    U.resize(dim);
    V.resize(dim);
    T.resize(dim);

    for (auto current_cell : _grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();
        auto loc = at(i, j);
        A.set_element(loc, loc, 2 * div);
        A.set_element(loc, at(i - 1, j), -inv_dx2);
        A.set_element(loc, at(i + 1, j), -inv_dx2);
        A.set_element(loc, at(i, j + 1), -inv_dy2);
        A.set_element(loc, at(i, j - 1), -inv_dy2);
        U.set_element(loc, loc, 1);
        V.set_element(loc, loc, 1);
        T.set_element(loc, loc, 1);
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

    auto zero_val = [at](SparseMatrix<Real> &M, std::vector<Real> &RHS, int i, int j) {
        M.set_element(at(i, j), at(i, j), 0);
        /*  M.set_element(at(i, j), at(i + 1, j), 0);
          M.set_element(at(i, j), at(i - 1, j), 0);
          M.set_element(at(i, j), at(i, j - 1), 0);
          M.set_element(at(i, j), at(i, j + 1), 0);*/
        RHS[at(i, j)] = 0;
    };

    for (auto &boundary : _boundaries) {
        for (auto &cell : *boundary->_cells) {
            int i = cell->i();
            int j = cell->j();
            int id = cell->id();
            auto loc = at(i, j);
            if (loc == 0) {
                int a = 4;
            }
            switch (boundary->get_type()) {
            case 0: {
                // Outlet

            } break;
            case 1: {
                // Inlet
                auto inlet_vel_u = static_cast<InletBoundary *>(boundary.get())->_inlet_U[id];
                auto inlet_vel_v = static_cast<InletBoundary *>(boundary.get())->_inlet_V[id];

                // Apply noslip for now
                if (cell->is_border(border_position::RIGHT)) {

                    // Velocity
                    V.set_element(loc, at(i + 1, j), -1);
                    U_RHS[loc] += 0.5 * inlet_vel_u;
                } else if (cell->is_border(border_position::LEFT)) {

                    // Velocity
                    V.set_element(loc, at(i - 1, j), -1);
                    U_RHS[loc] += 0.5 * inlet_vel_u;
                } else if (cell->is_border(border_position::TOP)) {

                    // Velocity
                    U.set_element(loc, at(i, j + 1), -1);
                    V_RHS[loc] += 0.5 * inlet_vel_v;
                } else if (cell->is_border(border_position::BOTTOM)) {

                    // Velocity
                    U.set_element(loc, at(i, j - 1), -1);
                    V_RHS[loc] += 0.5 * inlet_vel_v;
                }

            } break;
            case 2: {
                // NoSlip

                auto wall_vel_u = static_cast<NoSlipWallBoundary *>(boundary.get())->_wall_velocity[id];
                auto wall_vel_v = static_cast<NoSlipWallBoundary *>(boundary.get())->_wall_velocity[id];
                // Apply noslip for now
                if (cell->borders().size() == 1) {
                    if (cell->is_border(border_position::RIGHT)) {
                        // Velocity
                        V_RHS[loc] += wall_vel_v;
                        zero_val(U, U_RHS, i, j);
                        auto neighbor = cell->neighbour(border_position::RIGHT);
                        /*  if (is_inlet(neighbor->id())) {
                              V_RHS[loc] += static_cast<InletBoundary *>(_boundaries[BND_INLET_IDX].get())
                                                ->_inlet_V[neighbor->id()] / 2;
                          } else {
                              V.set_element(loc, at(i + 1, j), -1);
                          }*/
                        V.set_element(loc, at(i + 1, j), -1);
                    } else if (cell->is_border(border_position::LEFT)) {
                        // Velocity
                        V_RHS[loc] += wall_vel_v;
                        zero_val(U, U_RHS, i - 1, j);
                        auto neighbor = cell->neighbour(border_position::LEFT);
                        /*  if (is_inlet(neighbor->id())) {
                              V_RHS[loc] += static_cast<InletBoundary *>(_boundaries[BND_INLET_IDX].get())
                                                ->_inlet_V[neighbor->id()] / 2;
                          } else {
                              V.set_element(loc, at(i - 1, j), -1);
                          }*/
                        V.set_element(loc, at(i - 1, j), -1);
                    } else if (cell->is_border(border_position::TOP)) {
                        // Velocity
                        U_RHS[loc] += wall_vel_u;
                        zero_val(V, V_RHS, i, j);
                        auto neighbor = cell->neighbour(border_position::TOP);
                        /*  if (is_inlet(neighbor->id())) {
                              U_RHS[loc] += static_cast<InletBoundary *>(_boundaries[BND_INLET_IDX].get())
                                                ->_inlet_U[neighbor->id()] / 2;
                          } else {
                              U.set_element(loc, at(i, j + 1), -1);
                          }*/
                        U.set_element(loc, at(i, j + 1), -1);
                    } else if (cell->is_border(border_position::BOTTOM)) {
                        // Velocity
                        U_RHS[loc] += wall_vel_u;
                        zero_val(V, V_RHS, i, j - 1);
                        auto neighbor = cell->neighbour(border_position::BOTTOM);
                        /*  if (is_inlet(neighbor->id())) {
                              U_RHS[loc] += static_cast<InletBoundary *>(_boundaries[BND_INLET_IDX].get())
                                                ->_inlet_U[neighbor->id()] / 2;
                          } else {
                              U.set_element(loc, at(i, j - 1), -1);
                          }*/
                        U.set_element(loc, at(i, j - 1), -1);
                    }
                } else if (cell->borders().size() == 2) {
                    if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::TOP)) {

                        // Velocity
                        U_RHS[at(i - 1, j)] += wall_vel_u;
                        V_RHS[at(i, j - 1)] += wall_vel_v;
                        zero_val(U, U_RHS, i, j);
                        zero_val(V, V_RHS, i, j);
                        U.set_element(at(i - 1, j), at(i - 1, j + 1), -1);
                        V.set_element(at(i, j - 1), at(i + 1, j - 1), -1);
                    } else if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::BOTTOM)) {

                        // Velocity
                        U_RHS[at(i - 1, j)] += wall_vel_u;
                        V_RHS[at(i, j)] += wall_vel_v;
                        zero_val(U, U_RHS, i, j);
                        zero_val(V, V_RHS, i, j - 1);
                        U.set_element(at(i - 1, j), at(i - 1, j - 1), -1);
                        V.set_element(loc, at(i + 1, j), -1);
                    } else if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::TOP)) {

                        // Velocity
                        U_RHS[loc] += wall_vel_u;
                        V_RHS[at(i, j - 1)] += wall_vel_v;
                        zero_val(U, U_RHS, i - 1, j);
                        zero_val(V, V_RHS, i, j);
                        U.set_element(loc, at(i, j + 1), -1);
                        V.set_element(at(i, j - 1), at(i - 1, j - 1), -1);
                    } else if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::BOTTOM)) {

                        zero_val(U, U_RHS, i - 1, j);
                        zero_val(V, V_RHS, i, j - 1);
                        // Velocity
                        U_RHS[loc] += wall_vel_u;
                        V_RHS[loc] += wall_vel_v;
                        U.set_element(loc, at(i, j - 1), -1);
                        V.set_element(loc, at(i - 1, j), -1);
                    }
                } else {
                    U.set_element(loc, loc, 1);
                    V.set_element(loc, loc, 1);
                }
            } break;
            case 3: {
                // FreeSlip
                if (cell->borders().size() == 1) {
                    if (cell->is_border(border_position::RIGHT)) {
                        // Velocity
                        zero_val(U, U_RHS, i, j);
                    } else if (cell->is_border(border_position::LEFT)) {
                        zero_val(U, U_RHS, i - 1, j);
                    } else if (cell->is_border(border_position::TOP)) {
                        zero_val(V, V_RHS, i, j);
                    } else if (cell->is_border(border_position::BOTTOM)) {
                        zero_val(V, V_RHS, i, j - 1);
                    }
                }
            } break;
            default:
                break;
            }
        }
    }
    if (_field.calc_temp) {
        T_RHS.resize(dim);
        for (auto &boundary : _boundaries) {
            for (auto &cell : *boundary->_cells) {
                int i = cell->i();
                int j = cell->j();
                int id = cell->id();
                auto loc = at(i, j);
                switch (boundary->get_type()) {
                case 0: {
                    // Outlet
                    if (cell->is_border(border_position::RIGHT)) {
                        // T
                        T.set_element(loc, at(i + 1, j), 1);
                    } else if (cell->is_border(border_position::LEFT)) {
                        // T
                        T.set_element(loc, at(i - 1, j), 1);
                    } else if (cell->is_border(border_position::TOP)) {
                        // T
                        T.set_element(loc, at(i, j + 1), 1);
                    } else if (cell->is_border(border_position::BOTTOM)) {
                        // T
                        T.set_element(loc, at(i, j - 1), 1);
                    }
                } break;
                case 1: {
                    // Inlet
                    auto wt = static_cast<InletBoundary *>(boundary.get())->_wall_temperature[id];
                    T_RHS[loc] += wt;
                    if (cell->is_border(border_position::RIGHT)) {
                        // T
                        T.set_element(loc, at(i + 1, j), -1);
                    } else if (cell->is_border(border_position::LEFT)) {
                        // T
                        T.set_element(loc, at(i - 1, j), -1);
                    } else if (cell->is_border(border_position::TOP)) {
                        // T
                        T.set_element(loc, at(i, j + 1), -1);
                    } else if (cell->is_border(border_position::BOTTOM)) {
                        // T
                        T.set_element(loc, at(i, j - 1), -1);
                    }
                }

                break;
                case 2: {
                    // NoSlip
                    auto wt = static_cast<InletBoundary *>(boundary.get())->_wall_temperature[id];
                    if (cell->is_border(border_position::RIGHT)) {
                        // T
                        if (wt != -1) {
                            T_RHS[loc] += wt;
                            T.set_element(loc, at(i + 1, j), -1);
                        } else {
                            T.set_element(loc, at(i + 1, j), 1);
                        }
                    } else if (cell->is_border(border_position::LEFT)) {
                        // T
                        if (wt != -1) {
                            T_RHS[loc] += wt;
                            T.set_element(loc, at(i - 1, j), -1);
                        } else {
                            T.set_element(loc, at(i - 1, j), 1);
                        }
                    } else if (cell->is_border(border_position::TOP)) {
                        // T
                        if (wt != -1) {
                            T_RHS[loc] += wt;
                            T.set_element(loc, at(i, j + 1), -1);
                        } else {
                            T.set_element(loc, at(i, j + 1), 1);
                        }
                    } else if (cell->is_border(border_position::BOTTOM)) {
                        // T
                        if (wt != -1) {
                            T_RHS[loc] += wt;
                            T.set_element(loc, at(i, j - 1), -1);
                        } else {
                            T.set_element(loc, at(i, j - 1), 1);
                        }
                    }
                    if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::TOP)) {
                        // T
                        if (wt != -1) {
                            T_RHS[loc] += wt;
                            T.set_element(loc, at(i + 1, j), -0.5);
                            T.set_element(loc, at(i, j + 1), -0.5);
                        } else {
                        }
                    }
                    if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::BOTTOM)) {
                        // T
                        if (wt != -1) {
                            T_RHS[loc] += wt;
                            T.set_element(loc, at(i + 1, j), -0.5);
                            T.set_element(loc, at(i, j - 1), -0.5);
                        } else {
                        }
                    }
                    if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::TOP)) {
                        // T
                        if (wt != -1) {
                            T_RHS[loc] += wt;
                            T.set_element(loc, at(i - 1, j), -0.5);
                            T.set_element(loc, at(i, j + 1), -0.5);
                        } else {
                        }
                    }
                    if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::BOTTOM)) {
                        // T
                        if (wt != -1) {
                            T_RHS[loc] += wt;
                            T.set_element(loc, at(i - 1, j), -0.5);
                            T.set_element(loc, at(i, j - 1), -0.5);
                        } else {
                        }
                    }
                } break;
                case 3: {
                    // FreeSlip
                    auto wt = static_cast<InletBoundary *>(boundary.get())->_wall_temperature[id];
                    if (cell->is_border(border_position::RIGHT)) {
                        // T
                        if (wt != -1) {
                            T_RHS[loc] += wt;
                            T.set_element(loc, at(i + 1, j), -1);
                        } else {
                            T.set_element(loc, at(i + 1, j), 1);
                        }
                    } else if (cell->is_border(border_position::LEFT)) {
                        // T
                        if (wt != -1) {
                            T_RHS[loc] += wt;
                            T.set_element(loc, at(i - 1, j), -1);
                        } else {
                            T.set_element(loc, at(i - 1, j), 1);
                        }
                    } else if (cell->is_border(border_position::TOP)) {
                        // T
                        if (wt != -1) {
                            T_RHS[loc] += wt;
                            T.set_element(loc, at(i, j + 1), -1);
                        } else {
                            T.set_element(loc, at(i, j + 1), 1);
                        }
                    } else if (cell->is_border(border_position::BOTTOM)) {
                        if (wt != -1) {
                            T_RHS[loc] += wt;
                            T.set_element(loc, at(i, j - 1), -1);
                        } else {
                            T.set_element(loc, at(i, j - 1), 1);
                        }
                    }
                    if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::TOP)) {
                        // T
                        if (wt != -1) {
                            T_RHS[loc] += wt;
                            T.set_element(loc, at(i + 1, j), -0.5);
                            T.set_element(loc, at(i, j + 1), -0.5);
                        } else {
                        }
                    }
                    if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::BOTTOM)) {
                        if (wt != -1) {
                            T_RHS[loc] += wt;
                            T.set_element(loc, at(i + 1, j), -0.5);
                            T.set_element(loc, at(i, j - 1), -0.5);
                        } else {
                        }
                    }
                    if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::TOP)) {
                        if (wt != -1) {
                            T_RHS[loc] += wt;
                            T.set_element(loc, at(i - 1, j), -0.5);
                            T.set_element(loc, at(i, j + 1), -0.5);
                        } else {
                        }
                    }
                    if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::BOTTOM)) {
                        if (wt != -1) {
                            T_RHS[loc] += wt;
                            T.set_element(loc, at(i - 1, j), -0.5);
                            T.set_element(loc, at(i, j - 1), -0.5);
                        } else {
                        }
                    }
                } break;
                default:
                    break;
                }
            }
        }
    }
    U_fixed.construct_from_matrix(U);
    V_fixed.construct_from_matrix(V);
    T_fixed.construct_from_matrix(T);
}
