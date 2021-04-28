#include "Boundary.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::apply(Fields &field) {
    int imax = field.u_matrix().imax();
    int jmax = field.u_matrix().jmax();

    // For left and right boundary.
    std::vector<double> verticalBoundaryValue(jmax, 0);
    // For top and bottom boundary.
    std::vector<double> horizontalBoundaryValue(imax, 0);
    
    // Set boundary value for x-velocity matrix.
    // Vertical
    std::fill(verticalBoundaryValue.begin(), verticalBoundaryValue.end(), 0);
    field.u_matrix().set_col(verticalBoundaryValue, 0);
    field.u_matrix().set_col(verticalBoundaryValue, imax);
    // Horizontal
    horizontalBoundaryValue = field.u_matrix().get_row(1);
    std::for_each(horizontalBoundaryValue.begin(), horizontalBoundaryValue.end(), [](double& x){x = -1*x;});
    field.u_matrix().set_row(horizontalBoundaryValue, 0);
    horizontalBoundaryValue = field.u_matrix().get_row(jmax);
    std::for_each(horizontalBoundaryValue.begin(), horizontalBoundaryValue.end(), [](double& x){x = -1*x;});
    field.u_matrix().set_row(horizontalBoundaryValue, jmax+1);

    // Set boundary value for y-velocity matrix.
    // Horizontal
    std::fill(horizontalBoundaryValue.begin(), horizontalBoundaryValue.end(), 0);
    field.v_matrix().set_row(horizontalBoundaryValue, 0);
    field.v_matrix().set_row(horizontalBoundaryValue, jmax);
    // Vertical
    verticalBoundaryValue = field.v_matrix().get_col(1);
    std::for_each(verticalBoundaryValue.begin(), verticalBoundaryValue.end(), [](double& x){x = -1*x;});
    field.v_matrix().set_col(verticalBoundaryValue, 0);
    verticalBoundaryValue = field.v_matrix().get_col(imax);
    std::for_each(verticalBoundaryValue.begin(), verticalBoundaryValue.end(), [](double& x){x = -1*x;});
    field.v_matrix().set_col(verticalBoundaryValue, imax+1);

    // Set boundary value for pressure matrix.

    // Set boundary value for x-momentum flux matrix.

    // Set boundary value for y-momentum flux matrix.
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::apply(Fields &field) {}
