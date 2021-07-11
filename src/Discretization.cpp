#include "Discretization.hpp"

#include <assert.h>
#include <cmath>

std::vector<Real> Discretization::_dx; 
std::vector<Real> Discretization::_dy;
Real Discretization::_gamma = 0.0;

Discretization::Discretization(std::vector<Real> dx, std::vector<Real> dy, Real gamma) {
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
                     _dx[i];
    Real result_dc = _gamma *
                     (std::abs(U(i, j) + U(i + 1, j)) / 2 * (U(i, j) - U(i + 1, j)) / 2 -
                      std::abs(U(i - 1, j) + U(i, j)) / 2 * (U(i - 1, j) - U(i, j)) / 2) /
                     _dx[i];
    result += result_fd + result_dc;

    // dUV/dy
    result_fd = (interpolate(V, i, j, 1, 0) * interpolate(U, i, j, 0, 1) -
                 interpolate(V, i, j - 1, 1, 0) * interpolate(U, i, j, 0, -1)) /
                _dy[j];
    result_dc = _gamma *
                (std::abs(V(i, j) + V(i + 1, j)) / 2 * (U(i, j) - U(i, j + 1)) / 2 -
                 std::abs(V(i, j - 1) + V(i + 1, j - 1)) / 2 * (U(i, j - 1) - U(i, j)) / 2) /
                _dy[j];
    result += result_fd + result_dc;

    return result;
}

Real Discretization::convection_v(const Matrix<Real> &U, const Matrix<Real> &V, int i, int j) {
    // Convection term is split into finite differences and donor cell scheme, gamma is [0,1]
    Real result = 0.0;

    // dV^2/dy
    Real result_fd = (interpolate(V, i, j, 0, 1) * interpolate(V, i, j, 0, 1) -
                      interpolate(V, i, j, 0, -1) * interpolate(V, i, j, 0, -1)) /
                     _dy[j];
    Real result_dc = _gamma *
                     (std::abs(V(i, j) + V(i, j + 1)) / 2 * (V(i, j) - V(i, j + 1)) / 2 -
                      std::abs(V(i, j - 1) + V(i, j)) / 2 * (V(i, j - 1) - V(i, j)) / 2) /
                      _dy[j];
    result += result_fd + result_dc;

    // dUV/dx
    result_fd = (interpolate(U, i, j, 0, 1) * interpolate(V, i, j, 1, 0) -
                 interpolate(U, i - 1, j, 0, 1) * interpolate(V, i, j, -1, 0)) /
                _dx[i];
    result_dc = _gamma *
                (std::abs(U(i, j) + U(i, j + 1)) / 2 * (V(i, j) - V(i + 1, j)) / 2 -
                 std::abs(U(i - 1, j) + U(i - 1, j + 1)) / 2 * (V(i - 1, j) - V(i, j)) / 2) /
                _dx[i];
    result += result_fd + result_dc;

    return result;
}


Real Discretization::convection_uT(const Matrix<Real> &U, const Matrix<Real> &T, int i, int j) {
    Real result = 0.0;
    // d(uT)/dx
    Real result_fd = (U(i, j) * interpolate(T, i, j, 1, 0) - U(i - 1, j) * interpolate(T, i - 1, j, 1, 0)) / _dx[i];
    Real result_dc =
        _gamma * (std::abs(U(i, j)) * diff(T, i, j, 1, 0) - std::abs(U(i - 1, j)) * diff(T, i - 1, j, 1, 0)) / _dx[i];
    result += result_fd + result_dc;

    return result;
}

Real Discretization::convection_vT(const Matrix<Real> &V, const Matrix<Real> &T, int i, int j) {
    Real result = 0.0;
    // d(vT)/dy
    Real result_fd = (V(i, j) * interpolate(T, i, j, 0, 1) - V(i, j - 1) * interpolate(T, i, j - 1, 0, 1)) / _dy[j];
    Real result_dc =
        _gamma * (std::abs(V(i, j)) * diff(T, i, j, 0, 1) - std::abs(V(i, j - 1)) * diff(T, i, j - 1, 0, 1)) / _dy[j];
    result += result_fd + result_dc;

    return result;
}


Real Discretization::convection_uKEPS(const Matrix<Real> &U, const Matrix<Real> &T, int i, int j) {
    Real result = 0.0;
    constexpr int METHOD = 0;
    auto dx = _dx[i];
    auto dy = _dy[j];
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
    auto dx = _dx[i];
    auto dy = _dy[j];
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
    auto dx_left = (_dx[i-1] + _dx[i]) / 2;
    auto dx_right = (_dx[i] + _dx[i+1]) / 2;
    auto dy_bot = (_dy[j - 1] + _dy[j]) / 2;
    auto dy_top = (_dy[j] + _dy[j + 1]) / 2;

    auto dx = _dx[i];
    auto dy = _dy[j];

    auto i_diff = nu_i(i, j) * (P(i + 1, j) - P(i, j)) / dx_right - nu_i(i - 1, j) * (P(i, j) - P(i - 1, j)) / dx_left;
    auto j_diff = nu_j(i, j) * (P(i, j + 1) - P(i, j)) / dy_top - nu_j(i, j - 1) * (P(i, j) - P(i, j - 1)) / dy_bot;

    Real result = 1 / coeff * (i_diff / dx + j_diff / dy);

    return result;
}

Real Discretization::mean_strain_rate_squared(const Matrix<Real> &U, const Matrix<Real> &V, int i, int j) {
    Real result = 0;
    constexpr int METHOD = 0;
    if (METHOD == 0) {
        auto dx = _dx[i];
        auto dy = _dy[j];
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
        Real dudy = (interpolate(U, i - 1, j + 1, 1, 0) - interpolate(U, i - 1, j - 1, 1, 0)) / (2 * _dy[j]);
        Real dvdx = (interpolate(V, i + 1, j - 1, 0, 1) - interpolate(V, i - 1, j - 1, 0, 1)) / (2 * _dx[i]);
        Real dudx = (2 / 3) * (U(i, j) - U(i - 1, j)) / _dx[i];
        Real dvdy = (2 / 3) * (V(i, j) - V(i, j - 1)) / _dy[j];
        result = 2 * (dudy + dvdx + dudx + dvdy) * (dudy + dvdx + dudx + dvdy);
    }

    return result;
}

Real Discretization::diffusion(const Matrix<Real> &A, int i, int j) {
    Real dx_left = (_dx[i-1] + _dx[i]) / 2;
    Real dx_right = (_dx[i] + _dx[i+1]) / 2;

    Real dxx = ((A(i - 1, j) - A(i, j)) / dx_left - (A(i, j) - A(i + 1, j)) / dx_right) / _dx[i];

    Real dy_bot = (_dy[j - 1] + _dy[j]) / 2;
    Real dy_top = (_dy[j] + _dy[j + 1]) / 2;

    Real dyy = ((A(i, j - 1) - A(i, j)) / dy_bot - (A(i, j) - A(i, j + 1)) / dy_top) / _dy[j];

    Real result = dxx + dyy;
    return result;
}

Real Discretization::laplacian(const Matrix<Real> &P, int i, int j) {
    Real dx_left = (_dx[i-1] + _dx[i]) / 2;
    Real dx_right = (_dx[i] + _dx[i+1]) / 2;

    Real dxx = ((P(i - 1, j) - P(i, j)) / dx_left - (P(i, j) - P(i + 1, j)) / dx_right) / _dx[i];

    Real dy_bot = (_dy[j - 1] + _dy[j]) / 2;
    Real dy_top = (_dy[j] + _dy[j + 1]) / 2;

    Real dyy = ((P(i, j - 1) - P(i, j)) / dy_bot - (P(i, j) - P(i, j + 1)) / dy_top) / _dy[j];

    Real result = dxx + dyy;
    return result;
}

Real Discretization::sor_helper(const Matrix<Real> &P, int i, int j) {
    Real dx_left = (_dx[i-1] + _dx[i]) / 2;
    Real dx_right = (_dx[i] + _dx[i+1]) / 2;
    Real dy_bot = (_dy[j - 1] + _dy[j]) / 2;
    Real dy_top = (_dy[j] + _dy[j + 1]) / 2;

    Real result = P(i - 1, j) / (dx_left * _dx[i]) + P(i + 1, j) / (dx_right * _dx[i]) 
                + P(i, j - 1) / (dy_bot * _dy[j]) + P(i, j + 1) / (dy_top * _dy[j]);
       
    // Real result = (P(i + 1, j) + P(i - 1, j)) / (_dx[i] * _dx[i]) + (P(i, j + 1) + P(i, j - 1)) / (_dy[j] * _dy[j]);
    return result;
}

Real Discretization::interpolate(const Matrix<Real> &A, int i, int j, int i_offset, int j_offset) {
    // Interpolate only works between adjacent cells
    assert(i_offset < 2 && i_offset > -2);
    assert(j_offset < 2 && j_offset > -2);

    Real result = A(i, j);
    Real dx_full = (_dx[i] + _dx[i + i_offset]);
    Real dy_full = (_dy[j] + _dy[j + j_offset]);
    // std::cout << "dy[j]: " << _dy[j] << " dy[j+off]: " << _dy[j + j_offset] << " Alll: " << _dy[j + j_offset] / dy_full << std::endl;
    if (j_offset == 0) {
        result = _dx[i + i_offset] / dx_full * A(i, j) + _dx[i] / dx_full * A(i + i_offset, j);
    }
    else if (i_offset == 0) {
        result = _dy[j + j_offset] / dy_full * A(i, j) + _dy[j] / dy_full * A(i, j + j_offset);
        
    }
    else {
        Real result_x1 = _dx[i + i_offset] / dx_full * A(i, j) + _dx[i] / dx_full * A(i + i_offset, j);
        Real result_x2 = _dx[i + i_offset] / dx_full * A(i, j + j_offset) + _dx[i] / dx_full * A(i + i_offset, j + j_offset);
        result = result = _dy[j + j_offset] / dy_full * result_x1 + _dy[j] / dy_full * result_x2;
    }
    return result;
}

Real Discretization::diff(const Matrix<Real> &A, int i, int j, int i_offset, int j_offset) {
    Real result = (A(i, j) - A(i + i_offset, j + j_offset)) / 2;
    return result;
}
