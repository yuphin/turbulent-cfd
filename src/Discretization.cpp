#include "Discretization.hpp"

#include <assert.h>
#include <cmath>

double Discretization::_dx = 0.0;
double Discretization::_dy = 0.0;
double Discretization::_gamma = 0.0;

Discretization::Discretization(double dx, double dy, double gamma) {
    _dx = dx;
    _dy = dy;
    _gamma = gamma;
}

double Discretization::convection_u(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    // Convection term is split into finite differences and donor cell scheme, gamma is [0,1]
    double result = 0.0;

    // dU^2/dx
    double result_fd = (interpolate(U, i, j, 1, 0) * interpolate(U, i, j, 1, 0) -
                        interpolate(U, i, j, -1, 0) * interpolate(U, i, j, -1, 0)) /
                       _dx;
    double result_dc = _gamma *
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

double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    // Convection term is split into finite differences and donor cell scheme, gamma is [0,1]
    double result = 0.0;

    // dV^2/dy
    double result_fd = (interpolate(V, i, j, 0, 1) * interpolate(V, i, j, 0, 1) -
                        interpolate(V, i, j, 0, -1) * interpolate(V, i, j, 0, -1)) /
                       _dy;
    double result_dc = _gamma *
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

double Discretization::convection_uT(const Matrix<double> &U, const Matrix<double> &T, int i, int j) {
    double result = 0.0;
    // d(uT)/dx
    double result_fd = (U(i, j) * interpolate(T, i, j, 1, 0) - U(i - 1, j) * interpolate(T, i - 1, j, 1, 0)) / _dx;
    double result_dc =
        _gamma * (std::abs(U(i, j)) * diff(T, i, j, 1, 0) - std::abs(U(i - 1, j)) * diff(T, i - 1, j, 1, 0)) / _dx;
    result += result_fd + result_dc;

    return result;
}

double Discretization::convection_vT(const Matrix<double> &V, const Matrix<double> &T, int i, int j) {
    double result = 0.0;
    // d(vT)/dy
    double result_fd = (V(i, j) * interpolate(T, i, j, 0, 1) - V(i, j - 1) * interpolate(T, i, j - 1, 0, 1)) / _dy;
    double result_dc =
        _gamma * (std::abs(V(i, j)) * diff(T, i, j, 0, 1) - std::abs(V(i, j - 1)) * diff(T, i, j - 1, 0, 1)) / _dy;
    result += result_fd + result_dc;

    return result;
}

double Discretization::diffusion(const Matrix<double> &A, int i, int j) {
    // Same as laplacian?
    double result = (A(i + 1, j) - 2.0 * A(i, j) + A(i - 1, j)) / (_dx * _dx) +
                    (A(i, j + 1) - 2.0 * A(i, j) + A(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::laplacian(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) / (_dx * _dx) +
                    (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::sor_helper(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) + P(i - 1, j)) / (_dx * _dx) + (P(i, j + 1) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset) {
    double result = (A(i, j) + A(i + i_offset, j + j_offset)) / 2;
    return result;
}

double Discretization::diff(const Matrix<double> &A, int i, int j, int i_offset, int j_offset) {
    double result = (A(i, j) - A(i + i_offset, j + j_offset)) / 2;
    return result;
}