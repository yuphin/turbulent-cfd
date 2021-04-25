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

void Fields::calculate_fluxes(Grid &grid) {}

void Fields::calculate_rs(Grid &grid) {
    // Note: Index 0 and imax+1 are reserved for ghost cells
    for (int i = 1; i <= grid.imax(); i++) {
        for (int j = 1; j <= grid.jmax(); j++) {
            double f_diff = 1 / grid.dx() * (f(i, j) - f(i - 1, j));
            double g_diff = 1 / grid.dy() * (g(i, j) - g(i, j - 1));
            rs(i, j) = 1 / dt() * (f_diff + g_diff);
        }
    }
}

void Fields::calculate_velocities(Grid &grid) {
    for (int i = 0; i <= grid.imax(); i++) {
        for (int j = 0; j <= grid.jmax(); j++) {
            if (i <= grid.imax() - 1) {
                u(i, j) = f(i, j) - dt() / grid.dx() * (p(i + 1, j) - p(i, j));
            }
            if (j <= grid.jmax() - 1) {
                v(i, j) = g(i, j) - dt() / grid.dy() * (p(i, j + 1) - p(i, j));
            }
        }
    }
}

double Fields::calculate_dt(Grid &grid) { return _dt; }

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }

double Fields::dt() const { return _dt; }