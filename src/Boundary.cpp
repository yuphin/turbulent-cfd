#include "Boundary.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::apply(Fields &field) {
    // imax and jmax here include the ghost cells.
    int imax = field.u_matrix().imax()-2;
    int jmax = field.u_matrix().jmax()-2;


    // For left and right boundary.
    std::vector<double> verticalBoundaryValue(jmax+2, 0);
    // For bottom boundary.
    std::vector<double> horizontalBoundaryValue(imax+2, 0);
    
    // ===== Set boundary value for x-velocity matrix. =====
    // Vertical
    std::fill(verticalBoundaryValue.begin(), verticalBoundaryValue.end(), 0);
    // Left
    field.u_matrix().set_col(verticalBoundaryValue, 0);
    // Right
    field.u_matrix().set_col(verticalBoundaryValue, imax);
    
    // Horizontal
    // Bottom
    horizontalBoundaryValue = field.u_matrix().get_row(1);
    std::for_each(horizontalBoundaryValue.begin(), horizontalBoundaryValue.end(), [](double& x){x = -1.0*x;});
    field.u_matrix().set_row(horizontalBoundaryValue, 0);

    // ===== Set boundary value for y-velocity matrix. =====
    // Horizontal
    // Bottom
    std::fill(horizontalBoundaryValue.begin(), horizontalBoundaryValue.end(), 0);
    field.v_matrix().set_row(horizontalBoundaryValue, 0);
    
    // Vertical
    // Left
    verticalBoundaryValue = field.v_matrix().get_col(1);
    std::for_each(verticalBoundaryValue.begin(), verticalBoundaryValue.end(), [](double& x){x = -1.0*x;});
    field.v_matrix().set_col(verticalBoundaryValue, 0);
    // Right
    verticalBoundaryValue = field.v_matrix().get_col(imax);
    std::for_each(verticalBoundaryValue.begin(), verticalBoundaryValue.end(), [](double& x){x = -1.0*x;});
    field.v_matrix().set_col(verticalBoundaryValue, imax+1);

    // ===== Set boundary value for pressure matrix. =====
    // Horizontal
    // Bottom
    horizontalBoundaryValue = field.p_matrix().get_row(1);
    field.p_matrix().set_row(horizontalBoundaryValue, 0);
    
    // Vertical
    // Left
    verticalBoundaryValue = field.p_matrix().get_col(1);
    field.p_matrix().set_col(verticalBoundaryValue, 0);
    // Right
    verticalBoundaryValue = field.p_matrix().get_col(imax);
    field.p_matrix().set_col(verticalBoundaryValue, imax+1);

    // ===== Set boundary value for x-momentum flux matrix. =====
    // Left
    verticalBoundaryValue = field.u_matrix().get_col(0);
    field.f_matrix().set_col(verticalBoundaryValue, 0);
    // Right
    verticalBoundaryValue = field.u_matrix().get_col(imax);
    field.f_matrix().set_col(verticalBoundaryValue, imax);

    // ===== Set boundary value for y-momentum flux matrix. =====
    // Bottom
    horizontalBoundaryValue = field.v_matrix().get_row(0);
    field.g_matrix().set_row(horizontalBoundaryValue, 0);
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::apply(Fields &field) {
    // imax and jmax here include the ghost cells.
    int imax = field.u_matrix().imax()-2;
    int jmax = field.u_matrix().jmax()-2;

    // For top boundary.
    std::vector<double> horizontalBoundaryValue(imax+2, 0);
    
    // Set boundary value for x-velocity matrix.
    // Horizontal
    // Top
    horizontalBoundaryValue = field.u_matrix().get_row(jmax);
    double wall_vel = _wall_velocity[LidDrivenCavity::moving_wall_id];
    std::for_each(horizontalBoundaryValue.begin(), horizontalBoundaryValue.end(), [wall_vel](double& x){x = 2*wall_vel-x;});
    field.u_matrix().set_row(horizontalBoundaryValue, jmax+1);

    // Set boundary value for y-velocity matrix.
    // Horizontal
    // Top
    std::fill(horizontalBoundaryValue.begin(), horizontalBoundaryValue.end(), 0);
    field.v_matrix().set_row(horizontalBoundaryValue, jmax);

    // Set boundary value for pressure matrix.
    // Horizontal
    // Top
    horizontalBoundaryValue = field.p_matrix().get_row(jmax);
    field.p_matrix().set_row(horizontalBoundaryValue, jmax+1);

    // Set boundary value for y-momentum flux matrix.
    // Top
    horizontalBoundaryValue = field.v_matrix().get_row(jmax);
    field.g_matrix().set_row(horizontalBoundaryValue, jmax);
}
