#include "Boundary.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <mpi.h>

/////////// Constructors ///////////
OutletBoundary::OutletBoundary(std::vector<Cell *> *cells) : Boundary(cells) { _type = 0; }

InletBoundary::InletBoundary(std::vector<Cell *> *cells) : Boundary(cells) { _type = 1; }

InletBoundary::InletBoundary(std::vector<Cell *> *cells, std::unordered_map<int, Real> inlet_U,
                             std::unordered_map<int, Real> inlet_V, std::unordered_map<int, Real> inlet_T,
                             std::unordered_map<int, Real> inlet_K, std::unordered_map<int, Real> inlet_EPS, Real DP)
    : Boundary(cells), _inlet_U(inlet_U), _inlet_V(inlet_V), _inlet_T(inlet_T), _inlet_K(inlet_K),
      _inlet_EPS(inlet_EPS) {
    _type = 1;
}

NoSlipWallBoundary::NoSlipWallBoundary(std::vector<Cell *> *cells) : Boundary(cells) { _type = 2; }

NoSlipWallBoundary::NoSlipWallBoundary(std::vector<Cell *> *cells, std::unordered_map<int, Real> wall_velocity,
                                       std::unordered_map<int, Real> wall_temperature)
    : Boundary(cells), _wall_velocity(wall_velocity) {
    _wall_temperature = wall_temperature;
    _type = 2;
}

FreeSlipWallBoundary::FreeSlipWallBoundary(std::vector<Cell *> *cells) : Boundary(cells) { _type = 3; }

FreeSlipWallBoundary::FreeSlipWallBoundary(std::vector<Cell *> *cells, std::unordered_map<int, Real> wall_velocity,
                                           std::unordered_map<int, Real> wall_temperature)
    : Boundary(cells), _wall_velocity(wall_velocity) {
    _wall_temperature = wall_temperature;
    _type = 3;
}

/////////// Main Boundary ///////////
void Boundary::enforce_fg(Fields &field) {
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        if (cell->is_border(border_position::RIGHT)) {
            field.f(i, j) = field.u(i, j);
        }
        if (cell->is_border(border_position::LEFT)) {
            field.f(i - 1, j) = field.u(i - 1, j);
        }
        if (cell->is_border(border_position::TOP)) {
            field.g(i, j) = field.v(i, j);
        }
        if (cell->is_border(border_position::BOTTOM)) {
            field.g(i, j - 1) = field.v(i, j - 1);
        }
    }
}

void Boundary::enforce_p(Fields &field) {
    for (auto &cell : *_cells) {
        if (cell->borders().size() > 2) {
            std::cerr << "Forbidden cells!!" << std::endl;
            assert(false);
        } else if (cell->borders().size() == 2) {
            enforce_p_diagonal(field, cell);
        } else if (cell->borders().size() == 1) {
            enforce_p_main(field, cell);
        }
    }
}

void Boundary::enforce_p_main(Fields &field, Cell *cell) {
    int i = cell->i();
    int j = cell->j();
    if (cell->is_border(border_position::RIGHT)) {
        field.p(i, j) = field.p(i + 1, j);
    } else if (cell->is_border(border_position::LEFT)) {
        field.p(i, j) = field.p(i - 1, j);
    } else if (cell->is_border(border_position::TOP)) {
        field.p(i, j) = field.p(i, j + 1);
    } else if (cell->is_border(border_position::BOTTOM)) {
        field.p(i, j) = field.p(i, j - 1);
    }
}

void Boundary::enforce_p_diagonal(Fields &field, Cell *cell) {
    int i = cell->i();
    int j = cell->j();
    if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::TOP)) {
        field.p(i, j) = (field.p(i + 1, j) + field.p(i, j + 1)) / 2.0;
    }
    if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::BOTTOM)) {
        field.p(i, j) = (field.p(i + 1, j) + field.p(i, j - 1)) / 2.0;
    }
    if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::TOP)) {
        field.p(i, j) = (field.p(i - 1, j) + field.p(i, j + 1)) / 2.0;
    }
    if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::BOTTOM)) {
        field.p(i, j) = (field.p(i - 1, j) + field.p(i, j - 1)) / 2.0;
    }
}

void Boundary::enforce_t(Fields &field) {
    for (auto &cell : *_cells) {
        size_t num_borders = cell->borders().size();
        assert(num_borders <= 2);
        auto wt = _wall_temperature[cell->id()];
        if (num_borders == 1) {
            if (wt != -1) {
                enforce_t_drichlet_main(field, cell);
            } else {
                enforce_t_adiabatic_main(field, cell);
            }
        } else {
            if (wt != -1) {
                enforce_t_drichlet_diag(field, cell);
            } else {
                enforce_t_adiabatic_diag(field, cell);
            }
        }
    }
}

void Boundary::enforce_t_drichlet_main(Fields &field, Cell *cell) {
    int i = cell->i();
    int j = cell->j();
    auto wt = _wall_temperature[cell->id()];
    if (cell->is_border(border_position::RIGHT)) {
        field.t(i, j) = 2 * wt - field.t(i + 1, j);
    }
    if (cell->is_border(border_position::LEFT)) {
        field.t(i, j) = 2 * wt - field.t(i - 1, j);
    }
    if (cell->is_border(border_position::TOP)) {
        field.t(i, j) = 2 * wt - field.t(i, j + 1);
    }
    if (cell->is_border(border_position::BOTTOM)) {
        field.t(i, j) = 2 * wt - field.t(i, j - 1);
    }
}

void Boundary::enforce_t_drichlet_diag(Fields &field, Cell *cell) {
    int i = cell->i();
    int j = cell->j();
    auto wt = _wall_temperature[cell->id()];
    if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::TOP)) {
        field.t(i, j) = 2 * wt - (field.t(i + 1, j) + field.t(i, j + 1) / 2);
    } else if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::BOTTOM)) {
        field.t(i, j) = 2 * wt - (field.t(i + 1, j) + field.t(i, j - 1) / 2);
    }
    if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::TOP)) {
        field.t(i, j) = 2 * wt - (field.t(i - 1, j) + field.t(i, j + 1) / 2);
    }
    if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::BOTTOM)) {
        field.t(i, j) = 2 * wt - (field.t(i - 1, j) + field.t(i, j - 1) / 2);
    }
}

void Boundary::enforce_t_adiabatic_main(Fields &field, Cell *cell) {
    int i = cell->i();
    int j = cell->j();
    if (cell->is_border(border_position::RIGHT)) {
        field.t(i, j) = field.t(i + 1, j);
    }
    if (cell->is_border(border_position::LEFT)) {
        field.t(i, j) = field.t(i - 1, j);
    }
    if (cell->is_border(border_position::TOP)) {
        field.t(i, j) = field.t(i, j + 1);
    }
    if (cell->is_border(border_position::BOTTOM)) {
        field.t(i, j) = field.t(i, j - 1);
    }
}

void Boundary::enforce_t_adiabatic_diag(Fields &field, Cell *cell) {
    int i = cell->i();
    int j = cell->j();
    if (cell->is_border(border_position::RIGHT)) {
        field.t(i, j) = field.t(i + 1, j);
    }
    if (cell->is_border(border_position::LEFT)) {
        field.t(i, j) = field.t(i - 1, j);
    }
    if (cell->is_border(border_position::TOP)) {
        field.t(i, j) = field.t(i, j + 1);
    }
    if (cell->is_border(border_position::BOTTOM)) {
        field.t(i, j) = field.t(i, j - 1);
    }
}

void Boundary::enforce_nu_t(Fields &field) {
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        if (cell->is_border(border_position::RIGHT)) {
            field.k(i, j) = field.k(i + 1, j);
            field.eps(i, j) = field.eps(i + 1, j);
            field.nu_t(i, j) = field.nu_t(i + 1, j);
        }
        if (cell->is_border(border_position::LEFT)) {
            field.k(i, j) = field.k(i - 1, j);
            field.eps(i, j) = field.eps(i - 1, j);
            field.nu_t(i, j) = field.nu_t(i - 1, j);
        }
        if (cell->is_border(border_position::TOP)) {
            field.k(i, j) = field.k(i, j + 1);
            field.eps(i, j) = field.eps(i, j + 1);
            field.nu_t(i, j) = field.nu_t(i, j + 1);
        }
        if (cell->is_border(border_position::BOTTOM)) {
            field.k(i, j) = field.k(i, j - 1);
            field.eps(i, j) = field.eps(i, j - 1);
            field.nu_t(i, j) = field.nu_t(i, j - 1);
        }
    }
}

/////////// Outlet ///////////
void OutletBoundary::enforce_p(Fields &field) {
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        int id = cell->id();
        field.p(i, j) = field._PI;
    }
}

void OutletBoundary::enforce_t(Fields &field) {
    for (auto &cell : *_cells) {
        enforce_t_adiabatic_main(field, cell);
    }
}

/////////// Inlet ///////////
void InletBoundary::enforce_uv(Fields &field) {
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        int id = cell->id();

        if (cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = _inlet_U[id];
            field.v(i, j) = 2 * _inlet_V[id] - field.v(i + 1, j);
        } else if (cell->is_border(border_position::LEFT)) {
            field.u(i, j) = _inlet_U[id];
            field.v(i, j) = 2 * _inlet_V[id] - field.v(i - 1, j);
        } else if (cell->is_border(border_position::TOP)) {
            field.u(i, j) = 2 * _inlet_U[id] - field.u(i, j + 1);
            field.v(i, j) = _inlet_V[id];
        } else if (cell->is_border(border_position::BOTTOM)) {
            field.u(i, j) = 2 * _inlet_U[id] - field.u(i, j - 1);
            field.v(i, j) = _inlet_V[id];
        }
    }
}

void InletBoundary::enforce_p(Fields &field) {
    if (_inlet_DP == REAL_MAX) { // Inlet BC isn't specified, apply default
        Boundary::enforce_p(field);
        return;
    }
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        int id = cell->id();
        field.p(i, j) = _inlet_DP + field._PI;
    }
}

void InletBoundary::enforce_t(Fields &field) {
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        auto wt = _inlet_T[cell->id()];
        if (cell->is_border(border_position::RIGHT)) {
            field.t(i, j) = 2 * wt - field.t(i + 1, j);
        }
        if (cell->is_border(border_position::LEFT)) {
            field.t(i, j) = 2 * wt - field.t(i - 1, j);
        }
        if (cell->is_border(border_position::TOP)) {
            field.t(i, j) = 2 * wt - field.t(i, j + 1);
        }
        if (cell->is_border(border_position::BOTTOM)) {
            field.t(i, j) = 2 * wt - field.t(i, j - 1);
        }
    }
}

void InletBoundary::enforce_nu_t(Fields &field) {
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        auto wk = _inlet_K[cell->id()];
        auto weps = _inlet_EPS[cell->id()];
        if (cell->is_border(border_position::RIGHT)) {
            field.k(i, j) = 2 * wk - field.k(i + 1, j);
            field.eps(i, j) = 2 * weps - field.eps(i + 1, j);
            field.nu_t(i, j) = 0.09 * field.k(i, j) * field.k(i, j) / field.eps(i, j);
        }
        if (cell->is_border(border_position::LEFT)) {
            field.k(i, j) = 2 * wk - field.k(i + 1, j);
            field.eps(i, j) = 2 * weps - field.eps(i + 1, j);
            field.nu_t(i, j) = 0.09 * field.k(i, j) * field.k(i, j) / field.eps(i, j);
        }
        if (cell->is_border(border_position::TOP)) {
            field.k(i, j) = 2 * wk - field.k(i + 1, j);
            field.eps(i, j) = 2 * weps - field.eps(i + 1, j);
            field.nu_t(i, j) = 0.09 * field.k(i, j) * field.k(i, j) / field.eps(i, j);
        }
        if (cell->is_border(border_position::BOTTOM)) {
            field.k(i, j) = 2 * wk - field.k(i + 1, j);
            field.eps(i, j) = 2 * weps - field.eps(i + 1, j);
            field.nu_t(i, j) = 0.09 * field.k(i, j) * field.k(i, j) / field.eps(i, j);
        }
    }    
}

/////////// NoSlip Walls ///////////
void NoSlipWallBoundary::enforce_uv(Fields &field) {

    for (auto &cell : *_cells) {
        if (cell->borders().size() > 2) {
            std::cerr << "Forbidden cells!!" << std::endl;
            assert(false);
        } else if (cell->borders().size() == 2) {
            enforce_uv_diagonal(field, cell);
        } else if (cell->borders().size() == 1) {
            enforce_uv_main(field, cell);
        }
    }
}

void NoSlipWallBoundary::enforce_uv_main(Fields &field, Cell *cell) {
    int i = cell->i();
    int j = cell->j();
    int id = cell->id();

    if (cell->is_border(border_position::RIGHT)) {
        field.u(i, j) = 0;
        field.v(i, j) = 2 * _wall_velocity[id] - field.v(i + 1, j);
    } else if (cell->is_border(border_position::LEFT)) {
        field.u(i - 1, j) = 0;
        field.v(i, j) = 2 * _wall_velocity[id] - field.v(i - 1, j);
    } else if (cell->is_border(border_position::TOP)) {
        field.v(i, j) = 0;
        field.u(i, j) = 2 * _wall_velocity[id] - field.u(i, j + 1);
    } else if (cell->is_border(border_position::BOTTOM)) {
        field.v(i, j - 1) = 0;
        field.u(i, j) = 2 * _wall_velocity[id] - field.u(i, j - 1);
    }
}

void NoSlipWallBoundary::enforce_uv_diagonal(Fields &field, Cell *cell) {
    int i = cell->i();
    int j = cell->j();
    int id = cell->id();
    // This shouldn't happen since index 0 is a ghost cell
    assert(!(i == 0 && j == 0) && "Illegal index");

    if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::TOP)) {
        field.u(i, j) = 0;
        field.u(i - 1, j) = 2 * _wall_velocity[id] - field.u(i - 1, j + 1);
        field.v(i, j) = 0;
        field.v(i, j - 1) = 2 * _wall_velocity[id] - field.v(i + 1, j - 1);
    } else if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::BOTTOM)) {
        field.u(i, j) = 0;
        field.u(i - 1, j) = 2 * _wall_velocity[id] - field.u(i - 1, j - 1);
        field.v(i, j) = 2 * _wall_velocity[id] - field.v(i + 1, j);
        field.v(i, j - 1) = 0;
    } else if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::TOP)) {
        field.u(i, j) = 2 * _wall_velocity[id] - field.u(i, j + 1);
        field.u(i - 1, j) = 0;
        field.v(i, j) = 0;
        field.v(i, j - 1) = 2 * _wall_velocity[id] - field.v(i - 1, j - 1);
    } else if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::BOTTOM)) {
        field.u(i, j) = 2 * _wall_velocity[id] - field.u(i, j - 1);
        field.u(i - 1, j) = 0;
        field.v(i, j) = 2 * _wall_velocity[id] - field.v(i - 1, j);
        field.v(i, j - 1) = 0;
    }
}

/////////// FreeSlip Walls ///////////
void FreeSlipWallBoundary::enforce_uv(Fields &field) {
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        int id = cell->id();
        assert(cell->borders().size() <= 2);
        if (cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = 0;
        }
        if (cell->is_border(border_position::LEFT)) {
            field.u(i - 1, j) = 0;
        }
        if (cell->is_border(border_position::TOP)) {
            field.v(i, j) = 0;
        }
        if (cell->is_border(border_position::BOTTOM)) {
            field.v(i, j - 1) = 0;
        }
    }
}