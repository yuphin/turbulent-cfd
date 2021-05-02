#include "Boundary.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::enforce_uv(Fields &field, Grid &grid) {
    int imax = field.u_matrix().imax() - 2;
    int jmax = field.u_matrix().jmax() - 2;
    for (int i = 0; i <= imax; i++) {
        field.u(i, 0) = -field.u(i, 1);
        field.v(i, 0) = 0;
        field.v(i, jmax) = 0;
    }
    for (int j = 1; j <= jmax; j++) {
        field.u(0, j) = 0;
        field.u(imax, j) = 0;
        field.v(0, j) = -field.v(1, j);
        field.v(imax + 1, j) = -field.v(imax, j);
    }
}

void FixedWallBoundary::enforce_fg(Fields &field, Grid &grid) {
    int imax = grid.imax();
    int jmax = grid.jmax();
    for (int i = 0; i <= imax; i++) {
        field.g(i, 0) = field.v(i, 0);
        field.g(i, jmax) = field.v(i, jmax);
    }
    for (int j = 0; j <= jmax; j++) {
        field.f(0, j) = field.u(0, j);
        field.f(imax, j) = field.u(imax, j);
    }
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::enforce_uv(Fields &field, Grid &grid) {
    // imax and jmax here include the ghost cells.
    int imax = grid.imax();
    int jmax = grid.jmax();
    double wall_vel = _wall_velocity[LidDrivenCavity::moving_wall_id];
    for (int i = 0; i <= imax; i++) {
        field.u(i, jmax + 1) = 2 * wall_vel - field.u(i, jmax); // 1;
    }
}
