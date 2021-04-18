#include "Cell.hpp"

#include <iostream>

Cell::Cell(int i, int j, cell_type type) : _i(i), _j(j), _type(type) {}

Cell::Cell(int i, int j, cell_type type, int id) : _i(i), _j(j), _type(type), _id(id) {}

// borders Get and Set
bool Cell::is_border(border_position position) const { return _border.at(static_cast<int>(position)); };

const Cell *Cell::neighbour(border_position position) const { return _neighbours.at(static_cast<int>(position)); }

void Cell::set_neighbour(Cell *cell, border_position position) { _neighbours.at(static_cast<int>(position)) = cell; }

const std::vector<border_position> &Cell::borders() const { return _borders; }

void Cell::add_border(border_position border) {
    _border.at(static_cast<int>(border)) = true;
    _borders.push_back(border);
}

int Cell::i() const { return _i; }

int Cell::j() const { return _j; }

cell_type Cell::type() const { return _type; }

int Cell::wall_id() const { return _id; }