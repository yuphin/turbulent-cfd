#include "Boundary.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cassert>

OutletBoundary::OutletBoundary(std::vector<Cell *> *cells) : Boundary(cells) {}

InletBoundary::InletBoundary(std::vector<Cell *> *cells) : Boundary(cells) {}

InletBoundary::InletBoundary(std::vector<Cell *> *cells, std::unordered_map<int, double> inlet_U,
                                        std::unordered_map<int, double> inlet_V,
                                        std::unordered_map<int, double> inlet_T)
    : Boundary(cells), _inlet_U(inlet_U), _inlet_V(inlet_V), _inlet_T(inlet_T) {}

NoSlipWallBoundary::NoSlipWallBoundary(std::vector<Cell *> *cells) : Boundary(cells) {}

NoSlipWallBoundary::NoSlipWallBoundary(std::vector<Cell *> *cells, std::unordered_map<int, double> wall_velocity, 
                                        std::unordered_map<int, double> wall_temperature)
    : Boundary(cells), _wall_temperature(wall_temperature), _wall_velocity(wall_velocity) {}

FreeSlipWallBoundary::FreeSlipWallBoundary(std::vector<Cell *> *cells) : Boundary(cells) {}

FreeSlipWallBoundary::FreeSlipWallBoundary(std::vector<Cell *> *cells, std::unordered_map<int, double> wall_velocity,
                                       std::unordered_map<int, double> wall_temperature)
    : Boundary(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}


void OutletBoundary::enforce_uv(Fields &field) {assert(false);}
void InletBoundary::enforce_uv(Fields &field) {assert(false);}

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
        int i = cell->i();
        int j = cell->j();
        if (cell->is_border(border_position::RIGHT)) {
            field.p(i, j) = field.p(i + 1, j);
        }
        if (cell->is_border(border_position::LEFT)) {
            field.p(i, j) = field.p(i - 1, j);
        }
        if (cell->is_border(border_position::TOP)) {
            field.p(i, j) = field.p(i, j + 1);
        }
        if (cell->is_border(border_position::BOTTOM)) {
            field.p(i, j) = field.p(i, j - 1);
        }
    }
}

void NoSlipWallBoundary::enforce_uv(Fields &field) {
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        int id = cell->id();
        if (cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = 0;
            field.v(i, j) = 2 * _wall_velocity[id] - field.v(i + 1, j);
        }
        if (cell->is_border(border_position::LEFT)) {
            field.u(i - 1, j) = 0;
            field.v(i, j) = 2 * _wall_velocity[id] - field.v(i - 1, j);
        }
        if (cell->is_border(border_position::TOP)) {
            field.v(i, j) = 0;
            field.u(i, j) = 2 * _wall_velocity[id] - field.u(i, j + 1);
        }
        if (cell->is_border(border_position::BOTTOM)) {
            field.v(i, j - 1) = 0;
            // std::cout << _wall_velocity[id] << std::endl;
            field.u(i, j) = 2 * _wall_velocity[id] - field.u(i, j - 1);
        }
    }
}

void FreeSlipWallBoundary::enforce_uv(Fields &field) {
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        if (cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = 0;
            field.v(i, j) = 2 * _wall_velocity[LidDrivenCavity::moving_wall_id] - field.v(i + 1, j);
        }
        if (cell->is_border(border_position::LEFT)) {
            field.u(i - 1, j) = 0;
            field.v(i, j) = 2 * _wall_velocity[LidDrivenCavity::moving_wall_id] - field.v(i - 1, j);
        }
        if (cell->is_border(border_position::TOP)) {
            field.v(i, j) = 0;
            field.u(i, j) = 2 * _wall_velocity[LidDrivenCavity::moving_wall_id] - field.u(i, j + 1);
        }
        if (cell->is_border(border_position::BOTTOM)) {
            field.v(i, j - 1) = 0;
            field.u(i, j) = 2 * _wall_velocity[LidDrivenCavity::moving_wall_id] - field.u(i, j - 1);
        }
    }
}