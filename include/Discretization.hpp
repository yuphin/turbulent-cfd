#pragma once

#include "Datastructures.hpp"

/**
 * @brief Static discretization methods to modify the fields
 *
 */
class Discretization {
  public:
    Discretization() = default;

    /**
     * @brief Constructor to set the discretization parameters
     *
     * @param[in] cell size in x direction
     * @param[in] cell size in y direction
     * @param[in] upwinding coefficient
     */
    Discretization(Real dx, Real dy, Real gamma);

    /**
     * @brief Diffusion discretization in 2D using central differences
     *
     * @param[in] data to be discretized
     * @param[in] x index
     * @param[in] y index
     *
     */
    static Real diffusion(const Matrix<Real> &A, int i, int j);

    /**
     * @brief Convection in x direction using donor-cell scheme
     *
     * @param[in] x-velocity field
     * @param[in] y-velocity field
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static Real convection_u(const Matrix<Real> &U, const Matrix<Real> &V, int i, int j);

    /**
     * @brief Convection in y direction using donor-cell scheme
     *
     * @param[in] x-velocity field
     * @param[in] y-velocity field
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static Real convection_v(const Matrix<Real> &U, const Matrix<Real> &V, int i, int j);

    /**
     * @brief Convection of temperature in x direction using donor-cell scheme
     *
     * @param[in] x-velocity field
     * @param[in] Temperature
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static Real convection_uT(const Matrix<Real> &U, const Matrix<Real> &T, int i, int j);

    /**
     * @brief Convection of temperature in y direction using donor-cell scheme
     *
     * @param[in] y-velocity field
     * @param[in] Temperature
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static Real convection_vT(const Matrix<Real> &V, const Matrix<Real> &T, int i, int j);

    static Real convection_uKEPS(const Matrix<Real> &U, const Matrix<Real> &KEPS, int i, int j);
    static Real convection_vKEPS(const Matrix<Real> &V, const Matrix<Real> &KEPS, int i, int j);
    static Real laplacian_nu(const Matrix<Real> &P, Real nu, const Matrix<Real> &nu_i, const Matrix<Real> &nu_j, int i,
                             int j,
                             Real coeff = 1);
    static Real mean_strain_rate_squared(const Matrix<Real> &U, const Matrix<Real> &V, int i, int j);
    /**
     * @brief Laplacian term discretization using central difference
     *
     * @param[in] data to be discretized
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static Real laplacian(const Matrix<Real> &P, int i, int j);

    /**
     * @brief Terms of laplacian needed for SOR, i.e. excluding unknown value at
     * (i,j)
     *
     * @param[in] data to be discretized
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static Real sor_helper(const Matrix<Real> &P, int i, int j);

    /**
     * @brief Linear interpolation
     *
     * @param[in] data to be interpolated
     * @param[in] x index
     * @param[in] y index
     * @param[in] x offset
     * @param[in] y offset
     * @param[out] result
     *
     */
    static Real interpolate(const Matrix<Real> &A, int i, int j, int i_offset, int j_offset);

    /**
     * @brief Difference with offsets
     *
     * @param[in] data
     * @param[in] x index
     * @param[in] y index
     * @param[in] x offset
     * @param[in] y offset
     * @param[out] result
     */
    static Real diff(const Matrix<Real> &A, int i, int j, int i_offset, int j_offset);

    static Real _gamma;

  private:
    static Real _dx;
    static Real _dy;
   
};
