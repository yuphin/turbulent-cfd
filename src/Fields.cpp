#include "Fields.hpp"
#include "Communication.hpp"

#include <algorithm>
#include <iostream>
#include <math.h>

Fields::Fields(Real nu, Real dt, Real tau, int imax, int jmax, Real UI, Real VI, Real PI, Real TI, Real KI, Real EPSI, Real alpha,
               Real beta, Real gx, Real gy)
    : _nu(nu), _dt(dt), _tau(tau), _alpha(alpha), _beta(beta), _gx(gx), _gy(gy) {
    _U = Matrix<Real>(imax + 2, jmax + 2, UI);
    _V = Matrix<Real>(imax + 2, jmax + 2, VI);
    _P = Matrix<Real>(imax + 2, jmax + 2, PI);
    _T = Matrix<Real>(imax + 2, jmax + 2, TI);

    _F = Matrix<Real>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<Real>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<Real>(imax + 2, jmax + 2, 0.0);

    _NU_T = Matrix<Real>(imax + 2, jmax + 2, 0.0);
    _K = Matrix<Real>(imax + 2, jmax + 2, KI);
    _EPS = Matrix<Real>(imax + 2, jmax + 2, EPSI);

    _PI = PI;
    _TI = TI;
}

void Fields::calculate_nu_t(Grid &grid) {
    calculate_k(grid);
    calculate_epsilon(grid);    

    for (const auto &current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();

        nu_t(i, j) = 0.09 * k(i, j) * k(i, j) / eps(i, j);
    }
}

void Fields::calculate_k(Grid &grid) {
    auto K_OLD = _K;

    for (const auto &current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();

        k(i, j) = k(i, j) + _dt * ( (_nu + nu_t(i, j) / 1.0) * Discretization::laplacian(K_OLD, i, j) - 
                                        Discretization::convection_uT(_U, K_OLD, i, j) - Discretization::convection_vT(_V, K_OLD, i, j) -
                                        eps(i, j) + nu_t(i, j) * Discretization::mean_strain_rate_squared(_U, _V, i, j));
    }
}

void Fields::calculate_epsilon(Grid &grid) {
    auto EPS_OLD = _EPS;

    for (const auto &current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();

        eps(i, j) = eps(i, j) + _dt * ( (_nu + nu_t(i, j) / 1.3) * Discretization::laplacian(EPS_OLD, i, j) - 
                                        Discretization::convection_uT(_U, EPS_OLD, i, j) - Discretization::convection_vT(_V, EPS_OLD, i, j) -
                                        1.92 * eps(i, j) * eps(i, j) / k(i, j) + 
                                        1.44 * eps(i, j) / k(i, j) * nu_t(i, j) * Discretization::mean_strain_rate_squared(_U, _V, i, j));
    }
}

Real Fields::damp_f2(int i, int j) {
    return 1 - 0.3 * std::exp(- std::pow(k(i, j) * k(i, j) / (_nu * eps(i, j)), 2));
}

Real Fields::damp_fnu(int i, int j) {
    return std::exp(- 3.4 / (std::pow(k(i, j) * k(i, j) / (_nu * eps(i, j) * 50), 2)));
}

void Fields::calculate_fluxes(Grid &grid, bool calc_temp) {
    for (const auto &current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();
        f(i, j) = u(i, j) + dt() * ((_nu + Discretization::interpolate(nu_t_matrix(), i, j, 1, 0)) * 
                                    Discretization::diffusion(u_matrix(), i, j) -
                                    Discretization::convection_u(u_matrix(), v_matrix(), i, j));
        g(i, j) = v(i, j) + dt() * ((_nu +Discretization::interpolate(nu_t_matrix(), i, j, 0, 1)) * 
                                    Discretization::diffusion(v_matrix(), i, j) -
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
        Real f_diff = 1 / grid.dx() * (f(i, j) - f(i - 1, j));
        Real g_diff = 1 / grid.dy() * (g(i, j) - g(i, j - 1));
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
    auto T_old = _T;
    for (const auto &current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();
        t(i, j) =
            t(i, j) + _dt * (_alpha * Discretization::laplacian(T_old, i, j) -
                             Discretization::convection_uT(_U, T_old, i, j) - Discretization::convection_vT(_V, T_old, i, j));
    }
}

Real Fields::calculate_dt(Grid &grid, bool calc_temp) {

    Real dx2 = grid.dx() * grid.dx();
    Real dy2 = grid.dy() * grid.dy();

    // CFL conditions
    Real uMax = *std::max_element(_U.data(), _U.data() + _U.size());
    Real uMin = *std::min_element(_U.data(), _U.data() + _U.size());
    Real vMax = *std::max_element(_V.data(), _V.data() + _V.size());
    Real vMin = *std::min_element(_V.data(), _V.data() + _V.size());
    Real maxAbsU = fabs(uMax) > fabs(uMin) ? fabs(uMax) : fabs(uMin);
    Real maxAbsV = fabs(vMax) > fabs(vMin) ? fabs(vMax) : fabs(vMin);

    // Get the global maximums
    maxAbsU = Communication::reduce_all(maxAbsU, MPI_MAX);
    maxAbsV = Communication::reduce_all(maxAbsV, MPI_MAX);

    Real cond_2 = grid.dx() / maxAbsU;
    Real cond_3 = grid.dy() / maxAbsV;

    Real minimum = std::min(cond_2, cond_3);

    // viscosity limit
    if (_nu != 0.0) {
        Real cond_1 = 1.0 / (2.0 * _nu) * ((dx2 * dy2) / (dx2 + dy2));
        minimum = std::min(minimum, cond_1);
    }

    // thermal diffusitivity limit
    if (calc_temp) {
        Real cond_4 = 1 / (2 * _alpha * (1 / dx2 + 1 / dy2));
        minimum = std::min(minimum, cond_4);
    }
    _dt = _tau * minimum;
    return _dt;
}

Real &Fields::p(int i, int j) { return _P(i, j); }
Real &Fields::u(int i, int j) { return _U(i, j); }
Real &Fields::v(int i, int j) { return _V(i, j); }
Real &Fields::t(int i, int j) { return _T(i, j); }
Real &Fields::f(int i, int j) { return _F(i, j); }
Real &Fields::g(int i, int j) { return _G(i, j); }
Real &Fields::rs(int i, int j) { return _RS(i, j); }
Real &Fields::k(int i, int j) { return _K(i, j); }
Real &Fields::eps(int i, int j) { return _EPS(i, j); }
Real &Fields::nu_t(int i, int j) { return _NU_T(i, j); }

Matrix<Real> &Fields::p_matrix() { return _P; }
Matrix<Real> &Fields::u_matrix() { return _U; }
Matrix<Real> &Fields::v_matrix() { return _V; }
Matrix<Real> &Fields::t_matrix() { return _T; }
Matrix<Real> &Fields::f_matrix() { return _F; }
Matrix<Real> &Fields::g_matrix() { return _G; }
Matrix<Real> &Fields::k_matrix() { return _K; }
Matrix<Real> &Fields::eps_matrix() { return _EPS; }
Matrix<Real> &Fields::nu_t_matrix() { return _NU_T; }

Real Fields::dt() const { return _dt; }