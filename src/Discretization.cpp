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