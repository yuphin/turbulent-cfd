#include "Grid.hpp"
#include "Enums.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

Grid::Grid(std::string geom_name, Domain &domain) {

    _domain = domain;

    _cells = Matrix<Cell>(_domain.size_x + 2, _domain.size_y + 2);

    if (geom_name.compare("NONE")) {
        std::vector<std::vector<int>> geometry_data(_domain.domain_size_x + 2,
                                                    std::vector<int>(_domain.domain_size_y + 2, 0));
        parse_geometry_file(geom_name, geometry_data);
        assign_cell_types(geometry_data);
    } else {
        build_lid_driven_cavity();
    }
}

void Grid::build_lid_driven_cavity() {
    std::vector<std::vector<int>> geometry_data(_domain.domain_size_x + 2,
                                                std::vector<int>(_domain.domain_size_y + 2, 0));

    for (int i = 0; i < _domain.domain_size_x + 2; ++i) {
        for (int j = 0; j < _domain.domain_size_y + 2; ++j) {
            // Bottom, left and right walls: no-slip
            if (i == 0 || j == 0 || i == _domain.domain_size_x + 1) {
                geometry_data.at(i).at(j) = LidDrivenCavity::fixed_wall_id;
            }
            // Top wall: moving wall
            else if (j == _domain.domain_size_y + 1) {
                geometry_data.at(i).at(j) = LidDrivenCavity::moving_wall_id;
            }
        }
    }
    assign_cell_types(geometry_data);
}

void Grid::assign_cell_types(std::vector<std::vector<int>> &geometry_data) {

    int i = 0;
    int j = 0;

    for (int j_geom = _domain.jmin; j_geom < _domain.jmax; ++j_geom) {
        {
            i = 0;
        }
        for (int i_geom = _domain.imin; i_geom < _domain.imax; ++i_geom) {
            if (geometry_data.at(i_geom).at(j_geom) == 0) {
                _cells(i, j) = Cell(i, j, cell_type::FLUID);
                _fluid_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == LidDrivenCavity::moving_wall_id) {
                _cells(i, j) = Cell(i, j, cell_type::MOVING_WALL, geometry_data.at(i_geom).at(j_geom));
                _moving_wall_cells.push_back(&_cells(i, j));
            } else {
                if (i == 0 || j == 0 || i == _domain.size_x + 1 || j == _domain.size_y + 1) {
                    // Outer walls
                    _cells(i, j) = Cell(i, j, cell_type::FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                    _fixed_wall_cells.push_back(&_cells(i, j));
                }
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

void Grid::parse_geometry_file(std::string filedoc, std::vector<std::vector<int>> &geometry_data) {

    int numcols, numrows, depth;

    std::ifstream infile(filedoc);
    std::stringstream ss;
    std::string inputLine = "";

    // First line : version
    getline(infile, inputLine);
    if (inputLine.compare("P2") != 0) {
        std::cerr << "First line of the PGM file should be P2" << std::endl;
    }

    // Second line : comment
    getline(infile, inputLine);

    // Continue with a stringstream
    ss << infile.rdbuf();
    // Third line : size
    ss >> numrows >> numcols;
    // Fourth line : depth
    ss >> depth;

    //int array[numrows][numcols];

    // Following lines : data
    for (int col = numcols - 1; col > -1; --col) {
        for (int row = 0; row < numrows; ++row) {
            ss >> geometry_data[row][col];
        }
    }

    infile.close();
}

int Grid::imax() const { return _domain.size_x; }
int Grid::jmax() const { return _domain.size_y; }

int Grid::imaxb() const { return _domain.size_x + 2; }
int Grid::jmaxb() const { return _domain.size_y + 2; }

Cell Grid::cell(int i, int j) const { return _cells(i, j); }

double Grid::dx() const { return _domain.dx; }

double Grid::dy() const { return _domain.dy; }

const Domain &Grid::domain() const { return _domain; }

const std::vector<Cell *> &Grid::fluid_cells() const { return _fluid_cells; }

const std::vector<Cell *> &Grid::fixed_wall_cells() const { return _fixed_wall_cells; }

const std::vector<Cell *> &Grid::moving_wall_cells() const { return _moving_wall_cells; }
