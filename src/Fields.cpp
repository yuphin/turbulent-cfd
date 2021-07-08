#include "Fields.hpp"
#include "Communication.hpp"

#include <algorithm>
#include <assert.h>
#include <iostream>
#include <math.h>

Fields::Fields(Real nu, Real dt, Real tau, int imax, int jmax, Real UI, Real VI, Real PI, Real TI, Real KI, Real EPSI,
               Real alpha, Real beta, Real gx, Real gy)
    : _nu(nu), _dt(dt), _tau(tau), _alpha(alpha), _beta(beta), _gx(gx), _gy(gy) {
    _U = Matrix<Real>(imax + 2, jmax + 2, UI);
    _V = Matrix<Real>(imax + 2, jmax + 2, VI);
    _P = Matrix<Real>(imax + 2, jmax + 2, PI);
    _T = Matrix<Real>(imax + 2, jmax + 2, TI);

    _F = Matrix<Real>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<Real>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<Real>(imax + 2, jmax + 2, 0.0);

    _NU_T = Matrix<Real>(imax + 2, jmax + 2, 0.0);
    _NU_I = Matrix<Real>(imax + 2, jmax + 2, 0.0);
    _NU_J = Matrix<Real>(imax + 2, jmax + 2, 0.0);
    _K = Matrix<Real>(imax + 2, jmax + 2, KI);
    _EPS = Matrix<Real>(imax + 2, jmax + 2, EPSI);

    _PI = PI;
    _TI = TI;
    _UI = UI;
    _VI = VI;
    _KI = KI;
    _EPSI = EPSI;
    calc_temp = false;
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
Real &Fields::nu_i(int i, int j) { return _NU_I(i, j); }
Real &Fields::nu_j(int i, int j) { return _NU_J(i, j); }

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

void Fields::calculate_nu_t(Grid &grid) {
    for (const auto &current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();

        Real fnu_coeff = 1;
        auto kij = k(i, j);
        auto epsij = eps(i, j);
        nu_t(i, j) = fnu_coeff * 0.09 * kij * kij / epsij + _nu;
        assert(!isnan(nu_t(i, j)));
        assert(!isinf(nu_t(i, j)));
        assert(nu_t(i, j) > 0);
    }
    for (const auto &current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();
        Real fnu_coeff = 1;

        auto num_i = (k(i, j) + k(i + 1, j)) / 2;
        auto denom_i = (eps(i, j) + eps(i + 1, j)) / 2;
        auto num_j = (k(i, j) + k(i, j + 1)) / 2;
        auto denom_j = (eps(i, j) + eps(i, j + 1)) / 2;
        nu_i(i, j) = fnu_coeff * 0.09 * num_i * num_i / denom_i;
        nu_j(i, j) = fnu_coeff * 0.09 * num_j * num_j / denom_j;
    }
    calculate_k_and_epsilon(grid);
}

void Fields::calculate_k_and_epsilon(Grid &grid) {
    auto K_OLD = _K;
    auto EPS_OLD = _EPS;

    for (const auto &current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();
        Real f2_coeff = 1;
        auto nut = nu_t(i, j);
        auto kij = K_OLD(i, j);
        auto eij = EPS_OLD(i, j);
        auto k1_1 = Discretization::convection_uKEPS(_U, K_OLD, i, j);
        auto k1_2 = Discretization::convection_vKEPS(_V, K_OLD, i, j);
        auto e1_1 = Discretization::convection_uKEPS(_U, EPS_OLD, i, j);
        auto e1_2 = Discretization::convection_vKEPS(_V, EPS_OLD, i, j);

        auto k2 = Discretization::laplacian_nu(K_OLD, _nu, _NU_I, _NU_J, i, j);
        auto e2 = Discretization::laplacian_nu(EPS_OLD, _nu, _NU_I, _NU_J, i, j, 1.3);

        auto k3 = nut * Discretization::mean_strain_rate_squared(_U, _V, i, j);
        auto e3 = 1.44 * eij * k3 / kij;
        auto e4 = f2_coeff * 1.92 * eij * eij / kij;
        auto kij_new = kij + _dt * (-(k1_1 + k1_2) + k2 + k3 - eij);
        auto epsij_new = eij + _dt * (-(e1_1 + e1_2) + e2 + e3 - e4);
        k(i, j) = kij_new;
        eps(i, j) = epsij_new;
        assert(kij_new > 0);
        assert(epsij_new > 0);
    }
}

Real Fields::damp_f2(int i, int j) { return 1 - 0.3 * std::exp(-std::pow(k(i, j) * k(i, j) / (_nu * eps(i, j)), 2)); }

Real Fields::damp_fnu(int i, int j) {
    return std::exp(-3.4 / (std::pow(k(i, j) * k(i, j) / (_nu * eps(i, j) * 50), 2)));
}

void Fields::calculate_fluxes(Grid &grid, bool calc_temp, bool turbulent) {
    for (const auto &current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();
        Real nu_term1, nu_term2;
        if (turbulent) {
            nu_term1 = Discretization::interpolate(nu_t_matrix(), i, j, 1, 0);
            nu_term2 = Discretization::interpolate(nu_t_matrix(), i, j, 0, 1);
        } else {
            nu_term1 = nu_term2 = _nu;
        }
        f(i, j) = u(i, j) + dt() * (nu_term1 * Discretization::diffusion(u_matrix(), i, j) -
                                    Discretization::convection_u(u_matrix(), v_matrix(), i, j));
        g(i, j) = v(i, j) + dt() * (nu_term2 * Discretization::diffusion(v_matrix(), i, j) -
                                    Discretization::convection_v(u_matrix(), v_matrix(), i, j));
        if (calc_temp) {
            f(i, j) -= _beta * _dt / 2 * (t(i, j) + t(i + 1, j)) * _gx;
            g(i, j) -= _beta * _dt / 2 * (t(i, j) + t(i, j + 1)) * _gy;
        }
    }
}

// TODO add buoyancy term to the turbulent models
void Fields::calculate_temperatures(Grid &grid) {
    auto T_old = _T;
    for (const auto &current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();
        t(i, j) = t(i, j) + _dt * (_alpha * Discretization::laplacian(T_old, i, j) -
                                   Discretization::convection_uT(_U, T_old, i, j) -
                                   Discretization::convection_vT(_V, T_old, i, j));
    }
}

Real Fields::calculate_dt(Grid &grid, bool calc_temp, bool turbulent) {

    Real dx2 = grid.dx() * grid.dx();
    Real dy2 = grid.dy() * grid.dy();

    // CFL conditions
    Real uMax = *std::max_element(_U.data(), _U.data() + _U.size());
    Real uMin = *std::min_element(_U.data(), _U.data() + _U.size());
    Real vMax = *std::max_element(_V.data(), _V.data() + _V.size());
    Real vMin = *std::min_element(_V.data(), _V.data() + _V.size());
    Real maxAbsU = fabs(uMax) > fabs(uMin) ? fabs(uMax) : fabs(uMin);
    Real maxAbsV = fabs(vMax) > fabs(vMin) ? fabs(vMax) : fabs(vMin);
    Real nu_min = REAL_MAX;
    Real k_max;
    Real eps_max;
    if(turbulent){
        k_max = *std::max_element(_K.data(), _K.data() + _K.size());
        eps_max = *std::max_element(_EPS.data(), _EPS.data() + _EPS.size());
    }

    if (turbulent) {
        for (auto &cell : grid.fluid_cells()) {
            int i = cell->i();
            int j = cell->j();
            auto nut = nu_t(i, j);
            if (nut != 0 && nut < nu_min) {
                nu_min = nut;
            }
        }
    }
    nu_min = nu_min == REAL_MAX ? _nu : nu_min;

    // Get the global maximums
    maxAbsU = Communication::reduce_all(maxAbsU, MPI_MAX);
    maxAbsV = Communication::reduce_all(maxAbsV, MPI_MAX);

    Real cond_2 = grid.dx() / maxAbsU;
    Real cond_3 = grid.dy() / maxAbsV;

    Real minimum = std::min(cond_2, cond_3);

    // viscosity limit
    if (nu_min != 0.0) {
        Real cond_1 = 1.0 / (2.0 * _nu) * ((dx2 * dy2) / (dx2 + dy2));
        minimum = std::min(minimum, cond_1);
    }

    // thermal diffusitivity limit
    if (calc_temp) {
        Real cond_4 = 1 / (2 * _alpha * (1 / dx2 + 1 / dy2));
        minimum = std::min(minimum, cond_4);
    }
    if(turbulent){
        Real cond_5 = 1 / (2 * k_max * (1 / dx2 + 1 / dy2));
        Real cond_6 = 1 / (2 * eps_max * (1 / dx2 + 1 / dy2));
        minimum = std::min(minimum, cond_5);
        minimum = std::min(minimum, cond_6);
    }
    _dt = _tau * minimum;
    return _dt;
}