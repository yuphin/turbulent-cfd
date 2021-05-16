#include "Fields.hpp"

#include <algorithm>
#include <iostream>
#include <math.h>

Fields::Fields(double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double PI, double TI,
               double alpha, double beta, double gx, double gy)
    : _nu(nu), _dt(dt), _tau(tau), _alpha(alpha), _beta(beta), _gx(gx), _gy(gy) {
    _U = Matrix<double>(imax + 2, jmax + 2, UI);
    _V = Matrix<double>(imax + 2, jmax + 2, VI);
    _P = Matrix<double>(imax + 2, jmax + 2, PI);
    _T = Matrix<double>(imax + 2, jmax + 2, TI);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);
}

void Fields::calculate_fluxes(Grid &grid, bool calc_temp) {
    // Note: external forces e.g gravity not yet included
    for (const auto &current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();
        f(i, j) = u(i, j) + dt() * (_nu * Discretization::diffusion(u_matrix(), i, j) -
                                    Discretization::convection_u(u_matrix(), v_matrix(), i, j));
        g(i, j) = v(i, j) + dt() * (_nu * Discretization::diffusion(v_matrix(), i, j) -
                                    Discretization::convection_v(u_matrix(), v_matrix(), i, j));

        if (calc_temp) {
            f(i, j) -= _beta * _dt / 2 * (t(i, j) + t(i + 1, j)) * _gx;
            g(i, j) -= _beta * _dt / 2 * (t(i, j) + t(i, j + 1)) * _gy;
        }
    }
}

void Fields::calculate_rs(Grid &grid) {
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
        u(i, j) = f(i, j) - dt() / grid.dx() * (p(i + 1, j) - p(i, j));
        v(i, j) = g(i, j) - dt() / grid.dy() * (p(i, j + 1) - p(i, j));
    }
}

void Fields::calculate_temperatures(Grid &grid) {
    for (const auto &current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();
        t(i, j) = t(i, j) + _dt * (_alpha * Discretization::laplacian(_T, i, j) -
                             Discretization::convection_uT(_U, _T, i, j) - 
                             Discretization::convection_vT(_V, _T, i, j));
    }
}

double Fields::calculate_dt(Grid &grid, bool calc_temp) {
    double dx2 = grid.dx() * grid.dx();
    double dy2 = grid.dy() * grid.dy();

    // CFL conditions
    double uMax = *std::max_element(_U.data(), _U.data() + _U.size());
    double uMin = *std::min_element(_U.data(), _U.data() + _U.size());
    double vMax = *std::max_element(_V.data(), _V.data() + _V.size());
    double vMin = *std::min_element(_V.data(), _V.data() + _V.size());
    double maxAbsU = fabs(uMax) > fabs(uMin) ? fabs(uMax) : fabs(uMin);
    double maxAbsV = fabs(vMax) > fabs(vMin) ? fabs(vMax) : fabs(vMin);
    double cond_2 = grid.dx() / maxAbsU;
    double cond_3 = grid.dy() / maxAbsV;

    double minimum = std::min(cond_2, cond_3);

    // viscosity limit
    if (_nu != 0.0) {
        double cond_1 = 1.0 / (2.0 * _nu) * ((dx2 * dy2) / (dx2 + dy2));
        minimum = std::min(minimum, cond_1);
    }

    // thermal diffusitivity limit
    if (calc_temp) {
        double cond_4 = 1 / (2 * _alpha * (1 / dx2 + 1 / dy2));
        minimum = std::min(minimum, cond_4);
    }
    _dt = _tau * minimum;

    return _dt;
}

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::t(int i, int j) { return _T(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }
Matrix<double> &Fields::u_matrix() { return _U; }
Matrix<double> &Fields::v_matrix() { return _V; }
Matrix<double> &Fields::t_matrix() { return _T; }
Matrix<double> &Fields::f_matrix() { return _F; }
Matrix<double> &Fields::g_matrix() { return _G; }

double Fields::dt() const { return _dt; }