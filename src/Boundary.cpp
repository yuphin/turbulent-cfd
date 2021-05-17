#include "Boundary.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

OutletBoundary::OutletBoundary(std::vector<Cell *> *cells) : Boundary(cells) {}

InletBoundary::InletBoundary(std::vector<Cell *> *cells) : Boundary(cells) {}

InletBoundary::InletBoundary(std::vector<Cell *> *cells, std::unordered_map<int, double> inlet_U,
                                        std::unordered_map<int, double> inlet_V,
                                        std::unordered_map<int, double> inlet_T)
    : Boundary(cells), _inlet_U(inlet_U), _inlet_V(inlet_V), _inlet_T(inlet_T) {}

NoSlipWallBoundary::NoSlipWallBoundary(std::vector<Cell *> *cells) : Boundary(cells) {}

NoSlipWallBoundary::NoSlipWallBoundary(std::vector<Cell *> *cells, std::unordered_map<int, double> wall_velocity, 
                                        std::unordered_map<int, double> wall_temperature)
    : Boundary(cells),  _wall_velocity(wall_velocity) {
   _wall_temperature = wall_temperature;
}

FreeSlipWallBoundary::FreeSlipWallBoundary(std::vector<Cell *> *cells) : Boundary(cells) {}

FreeSlipWallBoundary::FreeSlipWallBoundary(std::vector<Cell *> *cells, std::unordered_map<int, double> wall_velocity,
                                       std::unordered_map<int, double> wall_temperature)
    : Boundary(cells), _wall_velocity(wall_velocity) {
    _wall_temperature = wall_temperature;
}

void Boundary::enforce_t_drichlet_main(Fields& field, Cell* cell) {
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
        field.t(i, j) = 3 * wt - (field.t(i + 1, j) + field.t(i, j + 1));
    } else if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::BOTTOM)) {
        field.t(i, j) = 3 * wt - (field.t(i + 1, j) + field.t(i, j - 1));
    }
    if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::TOP)) {
        field.t(i, j) = 3 * wt - (field.t(i - 1, j) + field.t(i, j + 1));
    }
    if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::BOTTOM)) {
        field.t(i, j) = 3 * wt - (field.t(i - 1, j) + field.t(i, j - 1));
    }
}

void Boundary::enforce_t_outflow_main(Fields &field, Cell *cell) {
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

void Boundary::enforce_t_outflow_diag(Fields &field, Cell *cell) {
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
void Boundary::enforce_t(Fields &field) {
	for (auto &cell : *_cells) {
		size_t num_borders = cell->borders().size();
		assert(num_borders <= 2);
        auto wt = _wall_temperature[cell->id()];
        if (num_borders == 1) {
			if (wt != -1) {
				enforce_t_drichlet_main(field, cell);
			} else {
                enforce_t_outflow_main(field, cell);
			}
        } else {
            if (wt != -1) {
                enforce_t_drichlet_diag(field, cell);
            } else {
                enforce_t_outflow_diag(field, cell);
            }
		}
	}
}

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
        if(cell->borders().size() >  2) {
            std::cerr << "Forbidden cells!!" << std::endl;
            assert(false);
        }
        else if(cell->borders().size() == 2) {
            enforce_p2(field, cell);
        }
        else if(cell->borders().size() == 1) {
            enforce_p1(field, cell);
        }           
    }
}

void Boundary::enforce_p1(Fields &field, Cell* cell) {
        int i = cell->i();
        int j = cell->j();
        if (cell->is_border(border_position::RIGHT)) {
            field.p(i, j) = field.p(i + 1, j);
        }
        else if (cell->is_border(border_position::LEFT)) {
            field.p(i, j) = field.p(i - 1, j);
        }
        else if (cell->is_border(border_position::TOP)) {
            field.p(i, j) = field.p(i, j + 1);
        }
        else if (cell->is_border(border_position::BOTTOM)) {
            field.p(i, j) = field.p(i, j - 1);
        }
}

void Boundary::enforce_p2(Fields &field, Cell* cell) {
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

void InletBoundary::enforce_uv(Fields &field) {
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        int id = cell->id();
        field.u(i, j) = _inlet_U[id];
        field.v(i, j) = _inlet_V[id];
        field.t(i, j) = _inlet_T[id];
    }
}

void NoSlipWallBoundary::enforce_uv(Fields &field) {
    for (auto &cell : *_cells) {
        if(cell->borders().size() >  2) {
            std::cerr << "Forbidden cells!!" << std::endl;
            assert(false);
        }
        else if(cell->borders().size() == 2) {
            enforce_uv2(field, cell);
        }
        else if(cell->borders().size() == 1) {
            enforce_uv1(field, cell);
        }           
    }
}

void NoSlipWallBoundary::enforce_uv2(Fields &field, Cell* cell) {
    int i = cell->i();
    int j = cell->j();
    int id = cell->id();
    
    if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::TOP)) {
        field.u(i, j) = 0;
        field.u(i - 1, j) = 2 * _wall_velocity[id] - field.u(i - 1, j + 1);
        field.v(i, j) = 0;
        field.v(i, j - 1) = 2 * _wall_velocity[id] - field.v(i + 1, j - 1);
    }
    else if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::BOTTOM)) {
        field.u(i, j) = 0;
        field.u(i - 1, j) = 2 * _wall_velocity[id] - field.u(i - 1, j - 1);
        field.v(i, j) = 2 * _wall_velocity[id] - field.v(i + 1, j);
        field.v(i, j - 1) = 0;
    }
    else if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::TOP)) {
        field.u(i, j) = 2 * _wall_velocity[id] - field.u(i, j + 1);
        field.u(i - 1, j) = 0;
        field.v(i, j) = 0;
        field.v(i, j - 1) = 2 * _wall_velocity[id] - field.v(i - 1, j - 1);
    }
    else if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::BOTTOM)) {
        field.u(i, j) = 2 * _wall_velocity[id] - field.u(i, j - 1);
        field.u(i - 1, j) = 0;
        field.v(i, j) = 2 * _wall_velocity[id] - field.v(i - 1, j);
        field.v(i, j - 1) = 0;
    }
}

void NoSlipWallBoundary::enforce_uv1(Fields &field, Cell* cell) {
    int i = cell->i();
    int j = cell->j();
    int id = cell->id();

    if (cell->is_border(border_position::RIGHT)) {
        field.u(i, j) = 0;
        field.v(i, j) = 2 * _wall_velocity[id] - field.v(i + 1, j);
        field.v(i, j - 1) = 2 * _wall_velocity[id] - field.v(i + 1, j - 1);
    }
    else if (cell->is_border(border_position::LEFT)) {
        field.u(i - 1, j) = 0;
        field.v(i, j) = 2 * _wall_velocity[id] - field.v(i - 1, j);
        field.v(i, j - 1) = 2 * _wall_velocity[id] - field.v(i - 1, j - 1);
    }
    else if (cell->is_border(border_position::TOP)) {
        field.v(i, j) = 0;
        field.u(i, j) = 2 * _wall_velocity[id] - field.u(i, j + 1);
        field.u(i - 1, j) = 2 * _wall_velocity[id] - field.u(i - 1, j + 1);
    }
    else if (cell->is_border(border_position::BOTTOM)) {
        field.v(i, j - 1) = 0;
        field.u(i, j) = 2 * _wall_velocity[id] - field.u(i, j - 1);
        field.u(i - 1, j) = 2 * _wall_velocity[id] - field.u(i - 1, j - 1);
    }
}

void FreeSlipWallBoundary::enforce_uv(Fields &field) {
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        int id = cell->id();
        if (cell->borders().size() > 2) {
            std::cerr << "Forbidden cells!!" << std::endl;
            assert(false);
        }
        if (cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = 0;
        } else if (cell->is_border(border_position::LEFT)) {
            field.u(i - 1, j) = 0;
        } else if (cell->is_border(border_position::TOP)) {
            field.v(i, j) = 0;
        } else if (cell->is_border(border_position::BOTTOM)) {
            field.v(i, j - 1) = 0;
        }
    }
}