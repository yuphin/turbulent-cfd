#include "Fields.hpp"

#include <algorithm>
#include <iostream>

Fields::Fields(double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double PI)
    : _nu(nu), _dt(dt), _tau(tau) {
    _U = Matrix<double>(imax + 2, jmax + 2, UI);
    _V = Matrix<double>(imax + 2, jmax + 2, VI);
    _P = Matrix<double>(imax + 2, jmax + 2, PI);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);
}

void Fields::calculate_fluxes(Grid &grid) {
    // Note: external forces e.g gravity not yet included
    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();
        f(i, j) = u(i, j) + dt() * (_nu * Discretization::diffusion(u_matrix(), i, j) 
                                    - Discretization::convection_u(u_matrix(), v_matrix(), i, j));
        g(i, j) = v(i, j) + dt() * (_nu * Discretization::diffusion(v_matrix(), i, j) 
                                    - Discretization::convection_v(u_matrix(), v_matrix(), i, j));
    }
}

void Fields::calculate_rs(Grid &grid) {
    // Note: Index 0 and imax+1 are reserved for ghost cells
    for (const auto &current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();
        double f_diff = 1 / grid.dx() * (f(i, j) - f(i - 1, j));
        double g_diff = 1 / grid.dy() * (g(i, j) - g(i, j - 1));
        rs(i, j) = 1 / dt() * (f_diff + g_diff);
    }
}

void Fields::calculate_velocities(Grid &grid) {
    for (const auto &current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();
        if (i <= grid.imax() - 1) {
            u(i, j) = f(i, j) - dt() / grid.dx() * (p(i + 1, j) - p(i, j));
        }
        if (j <= grid.jmax() - 1) {
            v(i, j) = g(i, j) - dt() / grid.dy() * (p(i, j + 1) - p(i, j));
        }
    }
}

double Fields::calculate_dt(Grid &grid) {
    double dx2 = grid.dx() * grid.dx(); 
    double dy2 = grid.dy() * grid.dy();
    double cond_1 = 1.0/(2.0*_nu) * ((dx2*dy2)/(dx2+dy2));

    double uMax = *std::max_element(_U.data(), _U.data()+_U.size());
    double vMax = *std::max_element(_V.data(), _V.data()+_V.size());
    double absUMax = std::abs(uMax);
    double absVMax = std::abs(vMax);
    double cond_2 = grid.dx() / absUMax;
    double cond_3 = grid.dy() / absVMax;

    double minimum = std::min(cond_1, cond_2);
    minimum = std::min(minimum, cond_3);
    _dt = _tau * minimum;

    return _dt;
}

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }
Matrix<double> &Fields::u_matrix() { return _U; }
Matrix<double> &Fields::v_matrix() { return _V; }


double Fields::dt() const { return _dt; }