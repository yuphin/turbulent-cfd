#include "Grid.hpp"
#include "Utilities.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

Grid::Grid(std::string geom_name, Domain &domain, std::vector<std::vector<int>> &geometry_data) {

    _domain = domain;

    _cells = Matrix<Cell>(_domain.size_x + 2, _domain.size_y + 2);
    assign_cell_types(geometry_data);
    preprocess_geometry();
}

void Grid::preprocess_geometry() {
    auto at = [this](int i, int j) { return &_cells._container[j * _domain.imax + i]; };
    auto get_idx = [this](int i, int j) { return j * _domain.imax + i; };
    for (auto cell : _fluid_cells) {
        int i = cell->i();
        int j = cell->j();
        cell_type ct;
        Real dist_x_p = 0;
        Real dist_y_p = 0;
        Real dist_x_n = 0;
        Real dist_y_n = 0;
        int idx_i_p = i;
        do {
            idx_i_p++;
            ct = at(idx_i_p, j)->type();
            // dist_x_p += at(idx_i_p, j)->dx;
            dist_x_p += dx();
        } while (ct == cell_type::FLUID);
        if (ct != cell_type::FREESLIP_WALL && ct != cell_type::NOSLIP_WALL) {
            dist_x_p = REAL_MAX;
        }
        int idx_i_n = i;
        do {
            idx_i_n--;
            ct = at(idx_i_n, j)->type();
            // dist_x_n += at(idx_i_n, j)->dx;
            dist_x_n += dx();
        } while (ct == cell_type::FLUID);
        if (ct != cell_type::FREESLIP_WALL && ct != cell_type::NOSLIP_WALL) {
            dist_x_n = REAL_MAX;
        }
        int idx_j_p = j;
        do {
            idx_j_p++;
            ct = at(i, idx_j_p)->type();
            // dist_y_p += at(i, idx_j_p)->dy;
            dist_y_p += dy();
        } while (ct == cell_type::FLUID);
        if (ct != cell_type::FREESLIP_WALL && ct != cell_type::NOSLIP_WALL) {
            dist_y_p = REAL_MAX;
        }
        int idx_j_n = j;
        do {
            idx_j_n--;
            ct = at(i, idx_j_n)->type();
            // dist_y_n += at(i, idx_j_n)->dy;
            dist_y_n += dy();
        } while (ct == cell_type::FLUID);
        if (ct != cell_type::FREESLIP_WALL && ct != cell_type::NOSLIP_WALL) {
            dist_y_n = REAL_MAX;
        }
        std::array<int, 4> idxs = {get_idx(idx_i_p, j), get_idx(idx_i_n, j), get_idx(i, idx_j_p), get_idx(i, idx_j_n)};
        std::array<Real, 4> dists = {dist_x_p, dist_x_n, dist_y_p, dist_y_n};
        Real *min_xy = std::min_element(dists.data(), dists.data() + 4);
        int min_idx = min_xy - dists.data();
        cell->closest_dist = *min_xy;
        cell->closest_wall_idx = idxs[min_idx];
    }
}

void Grid::assign_cell_types(std::vector<std::vector<int>> &geometry_data) {

    int i = 0;
    int j = 0;

    for (int j_geom = _domain.jmin; j_geom < _domain.jmax; ++j_geom) {
        { i = 0; }
        for (int i_geom = _domain.imin; i_geom < _domain.imax; ++i_geom) {
            if (geometry_data.at(i_geom).at(j_geom) == 0) {
                if (i_geom == _domain.imin || i_geom == _domain.imax - 1 
                    || j_geom == _domain.jmin || j_geom == _domain.jmax -1) {
                    ++i;
                    continue;
                }
                _cells(i, j) = Cell(i, j, cell_type::FLUID);
                _fluid_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == 1) {
                _cells(i, j) = Cell(i, j, cell_type::OUTLET);
                _outlet_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) >= 2 && geometry_data.at(i_geom).at(j_geom) <= 9) {
                _cells(i, j) = Cell(i, j, cell_type::INLET, geometry_data.at(i_geom).at(j_geom));
                _inlet_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) >= 10 && geometry_data.at(i_geom).at(j_geom) <= 19) {
                _cells(i, j) = Cell(i, j, cell_type::NOSLIP_WALL, geometry_data.at(i_geom).at(j_geom));
                _noslip_wall_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) >= 20) {
                _cells(i, j) = Cell(i, j, cell_type::FREESLIP_WALL, geometry_data.at(i_geom).at(j_geom));
                _freeslip_wall_cells.push_back(&_cells(i, j));
            }

            ++i;
        }
        ++j;
    }

    // Corner cell neighbour assigment
    // Bottom-Left Corner
    i = 0;
    j = 0;
    _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
    _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
    if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
        _cells(i, j).add_border(border_position::TOP);
    }
    if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
        _cells(i, j).add_border(border_position::RIGHT);
    }
    // Top-Left Corner
    i = 0;
    j = _domain.size_y + 1;
    _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
    _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
    if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
        _cells(i, j).add_border(border_position::BOTTOM);
    }
    if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
        _cells(i, j).add_border(border_position::RIGHT);
    }

    // Top-Right Corner
    i = _domain.size_x + 1;
    j = Grid::_domain.size_y + 1;
    _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
    _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
    if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
        _cells(i, j).add_border(border_position::BOTTOM);
    }
    if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
        _cells(i, j).add_border(border_position::LEFT);
    }

    // Bottom-Right Corner
    i = Grid::_domain.size_x + 1;
    j = 0;
    _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
    _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
    if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
        _cells(i, j).add_border(border_position::TOP);
    }
    if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
        _cells(i, j).add_border(border_position::LEFT);
    }
    // Bottom cells
    j = 0;
    for (int i = 1; i < _domain.size_x + 1; ++i) {
        _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
        _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
        _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::RIGHT);
        }
        if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::LEFT);
        }
        if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::TOP);
        }
    }

    // Top Cells
    j = Grid::_domain.size_y + 1;

    for (int i = 1; i < _domain.size_x + 1; ++i) {
        _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
        _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
        _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::RIGHT);
        }
        if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::LEFT);
        }
        if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::BOTTOM);
        }
    }

    // Left Cells
    i = 0;
    for (int j = 1; j < _domain.size_y + 1; ++j) {
        _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
        _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
        _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::RIGHT);
        }
        if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::BOTTOM);
        }
        if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::TOP);
        }
    }
    // Right Cells
    i = Grid::_domain.size_x + 1;
    for (int j = 1; j < _domain.size_y + 1; ++j) {
        _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
        _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
        _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
        if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::LEFT);
        }
        if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::BOTTOM);
        }
        if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::TOP);
        }
    }

    // Inner cells
    for (int i = 1; i < _domain.size_x + 1; ++i) {
        for (int j = 1; j < _domain.size_y + 1; ++j) {
            _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
            _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
            _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
            _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);

            if (_cells(i, j).type() != cell_type::FLUID) {
                if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
                    _cells(i, j).add_border(border_position::LEFT);
                }
                if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
                    _cells(i, j).add_border(border_position::RIGHT);
                }
                if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
                    _cells(i, j).add_border(border_position::BOTTOM);
                }
                if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
                    _cells(i, j).add_border(border_position::TOP);
                }
            }
        }
    }
}

int Grid::imax() const { return _domain.size_x; }
int Grid::jmax() const { return _domain.size_y; }

int Grid::imaxb() const { return _domain.size_x + 2; }
int Grid::jmaxb() const { return _domain.size_y + 2; }

Cell Grid::cell(int i, int j) const { return _cells(i, j); }

Real Grid::dx() const { return _domain.dx; }

Real Grid::dy() const { return _domain.dy; }

const Domain &Grid::domain() const { return _domain; }

const std::vector<Cell *> &Grid::fluid_cells() const { return _fluid_cells; }

std::vector<Cell *> &Grid::inlet_cells() { return _inlet_cells; }

std::vector<Cell *> &Grid::outlet_cells() { return _outlet_cells; }

std::vector<Cell *> &Grid::freeslip_wall_cells() { return _freeslip_wall_cells; }

std::vector<Cell *> &Grid::noslip_wall_cells() { return _noslip_wall_cells; }
