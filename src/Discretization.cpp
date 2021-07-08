#include "Discretization.hpp"

#include <assert.h>
#include <cmath>

Real Discretization::_dx = 0.0;
Real Discretization::_dy = 0.0;
Real Discretization::_gamma = 0.0;

Discretization::Discretization(Real dx, Real dy, Real gamma) {
    _dx = dx;
    _dy = dy;
    _gamma = gamma;
}

Real Discretization::convection_u(const Matrix<Real> &U, const Matrix<Real> &V, int i, int j) {
    // Convection term is split into finite differences and donor cell scheme, gamma is [0,1]
    Real result = 0.0;

    // dU^2/dx
    Real result_fd = (interpolate(U, i, j, 1, 0) * interpolate(U, i, j, 1, 0) -
                      interpolate(U, i, j, -1, 0) * interpolate(U, i, j, -1, 0)) /
                     _dx;
    Real result_dc = _gamma *
                     (std::abs(U(i, j) + U(i + 1, j)) / 2 * (U(i, j) - U(i + 1, j)) / 2 -
                      std::abs(U(i - 1, j) + U(i, j)) / 2 * (U(i - 1, j) - U(i, j)) / 2) /
                     _dx;
    result += result_fd + result_dc;

    // dUV/dy
    result_fd = (interpolate(V, i, j, 1, 0) * interpolate(U, i, j, 0, 1) -
                 interpolate(V, i, j - 1, 1, 0) * interpolate(U, i, j, 0, -1)) /
                _dy;
    result_dc = _gamma *
                (std::abs(V(i, j) + V(i + 1, j)) / 2 * (U(i, j) - U(i, j + 1)) / 2 -
                 std::abs(V(i, j - 1) + V(i + 1, j - 1)) / 2 * (U(i, j - 1) - U(i, j)) / 2) /
                _dy;
    result += result_fd + result_dc;

    return result;
}

Real Discretization::convection_v(const Matrix<Real> &U, const Matrix<Real> &V, int i, int j) {
    // Convection term is split into finite differences and donor cell scheme, gamma is [0,1]
    Real result = 0.0;

    // dV^2/dy
    Real result_fd = (interpolate(V, i, j, 0, 1) * interpolate(V, i, j, 0, 1) -
                      interpolate(V, i, j, 0, -1) * interpolate(V, i, j, 0, -1)) /
                     _dy;
    Real result_dc = _gamma *
                     (std::abs(V(i, j) + V(i, j + 1)) / 2 * (V(i, j) - V(i, j + 1)) / 2 -
                      std::abs(V(i, j - 1) + V(i, j)) / 2 * (V(i, j - 1) - V(i, j)) / 2) /
                     _dy;
    result += result_fd + result_dc;

    // dUV/dx
    result_fd = (interpolate(U, i, j, 0, 1) * interpolate(V, i, j, 1, 0) -
                 interpolate(U, i - 1, j, 0, 1) * interpolate(V, i, j, -1, 0)) /
                _dx;
    result_dc = _gamma *
                (std::abs(U(i, j) + U(i, j + 1)) / 2 * (V(i, j) - V(i + 1, j)) / 2 -
                 std::abs(U(i - 1, j) + U(i - 1, j + 1)) / 2 * (V(i - 1, j) - V(i, j)) / 2) /
                _dx;
    result += result_fd + result_dc;

    return result;
}

Real Discretization::convection_uT(const Matrix<Real> &U, const Matrix<Real> &T, int i, int j) {
    Real result = 0.0;
    // d(uT)/dx
    Real result_fd = (U(i, j) * interpolate(T, i, j, 1, 0) - U(i - 1, j) * interpolate(T, i - 1, j, 1, 0)) / _dx;
    Real result_dc =
        _gamma * (std::abs(U(i, j)) * diff(T, i, j, 1, 0) - std::abs(U(i - 1, j)) * diff(T, i - 1, j, 1, 0)) / _dx;
    result += result_fd + result_dc;

    return result;
}

Real Discretization::convection_vT(const Matrix<Real> &V, const Matrix<Real> &T, int i, int j) {
    Real result = 0.0;
    // d(vT)/dy
    Real result_fd = (V(i, j) * interpolate(T, i, j, 0, 1) - V(i, j - 1) * interpolate(T, i, j - 1, 0, 1)) / _dy;
    Real result_dc =
        _gamma * (std::abs(V(i, j)) * diff(T, i, j, 0, 1) - std::abs(V(i, j - 1)) * diff(T, i, j - 1, 0, 1)) / _dy;
    result += result_fd + result_dc;

    return result;
}

Real Discretization::convection_uKEPS(const Matrix<Real> &U, const Matrix<Real> &T, int i, int j) {
    Real result = 0.0;
    constexpr int METHOD = 0;
    auto dx = _dx;
    if (METHOD == 0) {
        Real result_fd = (U(i, j) * interpolate(T, i, j, 1, 0) - U(i - 1, j) * interpolate(T, i - 1, j, 1, 0)) / dx;
        Real result_dc =
            (std::abs(U(i, j)) * diff(T, i, j, 1, 0) - std::abs(U(i - 1, j)) * diff(T, i - 1, j, 1, 0)) / dx;
        result += result_fd + result_dc;
    } else {
        auto idnr = [&](int i, int j) -> Real { return U(i, j) <= 0; };
        result = (U(i, j) * T(i + idnr(i, j), j) - U(i - 1, j) * T(i - 1 + idnr(i - 1, j), j)) / dx;
    }
    return result;
}

Real Discretization::convection_vKEPS(const Matrix<Real> &V, const Matrix<Real> &T, int i, int j) {
    Real result = 0.0;
    constexpr int METHOD = 0;
    auto dy = _dy;
    if (METHOD == 0) {
        Real result_fd = (V(i, j) * interpolate(T, i, j, 0, 1) - V(i, j - 1) * interpolate(T, i, j - 1, 0, 1)) / dy;
        Real result_dc =
            1 * (std::abs(V(i, j)) * diff(T, i, j, 0, 1) - std::abs(V(i, j - 1)) * diff(T, i, j - 1, 0, 1)) / dy;
        result += result_fd + result_dc;
    } else {
        auto jdnr = [&](int i, int j) -> Real { return V(i, j) <= 0; };
        result = (V(i, j) * T(i, jdnr(i, j) + j) - V(i, j - 1) * T(i, j - 1 + jdnr(i, j - 1))) / dy;
    }
    return result;
}

Real Discretization::laplacian_nu(const Matrix<Real> &P, Real nu, const Matrix<Real> &nu_i, const Matrix<Real> &nu_j,
                                  int i, int j, Real coeff) {
    auto dx = _dx;
    auto dy = _dy;
    auto invdx = 1 / dx;
    auto invdx2 = invdx * invdx;
    auto invdy = 1 / dy;
    auto invdy2 = invdy * invdy;
    auto i_diff = (nu + nu_i(i, j)) * (P(i + 1, j) - P(i, j)) - (nu + nu_i(i - 1, j)) * (P(i, j) - P(i - 1, j));
    auto j_diff = (nu + nu_j(i, j)) * (P(i, j + 1) - P(i, j)) - (nu + nu_j(i, j - 1)) * (P(i, j) - P(i, j - 1));
    Real result = 1 / coeff * (i_diff * invdx2 + j_diff * invdy2);
    return result;
}

Real Discretization::mean_strain_rate_squared(const Matrix<Real> &U, const Matrix<Real> &V, int i, int j) {
    Real result = 0;
    constexpr int METHOD = 0;
    if (METHOD == 0) {
        auto dx = _dx;
        auto dy = _dy;
        auto invdx = 1 / dx;
        auto invdx2 = invdx * invdx;
        auto invdy = 1 / dy;
        auto invdy2 = invdy * invdy;
        Real u_diff = U(i, j) - U(i - 1, j);
        Real v_diff = V(i, j) - V(i, j - 1);
        auto shear_1 = (U(i, j + 1) + U(i - 1, j + 1) - U(i, j - 1) - U(i - 1, j - 1)) * (0.25 * invdy);
        auto shear_2 = (V(i + 1, j) + V(i + 1, j - 1) - V(i - 1, j) - V(i - 1, j - 1)) * (0.25 * invdx);
        auto shear = shear_1 + shear_2;
        result = (u_diff * u_diff) * invdx2 + v_diff * v_diff * invdy2 + shear * shear;
    } else {
        Real dudy = (interpolate(U, i - 1, j + 1, 1, 0) - interpolate(U, i - 1, j - 1, 1, 0)) / (2 * _dy);
        Real dvdx = (interpolate(V, i + 1, j - 1, 0, 1) - interpolate(V, i - 1, j - 1, 0, 1)) / (2 * _dx);
        Real dudx = (2 / 3) * (U(i, j) - U(i - 1, j)) / _dx;
        Real dvdy = (2 / 3) * (V(i, j) - V(i, j - 1)) / _dx;
        result = 2 * (dudy + dvdx + dudx + dvdy) * (dudy + dvdx + dudx + dvdy);
    }

    return result;
}

Real Discretization::diffusion(const Matrix<Real> &A, int i, int j) {
    // Same as laplacian?
    Real result = (A(i + 1, j) - 2.0 * A(i, j) + A(i - 1, j)) / (_dx * _dx) +
                  (A(i, j + 1) - 2.0 * A(i, j) + A(i, j - 1)) / (_dy * _dy);
    return result;
}

Real Discretization::laplacian(const Matrix<Real> &P, int i, int j) {
    Real result = (P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) / (_dx * _dx) +
                  (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

Real Discretization::sor_helper(const Matrix<Real> &P, int i, int j) {
    Real result = (P(i + 1, j) + P(i - 1, j)) / (_dx * _dx) + (P(i, j + 1) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

Real Discretization::interpolate(const Matrix<Real> &A, int i, int j, int i_offset, int j_offset) {
    Real result = (A(i, j) + A(i + i_offset, j + j_offset)) / 2;
    return result;
}

Real Discretization::diff(const Matrix<Real> &A, int i, int j, int i_offset, int j_offset) {
    Real result = (A(i, j) - A(i + i_offset, j + j_offset)) / 2;
    return result;
}