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
     * @param[in] dx cell size in x direction
     * @param[in] dy cell size in y direction
     * @param[in] gamma upwinding coefficient
     */
    Discretization(Real dx, Real dy, Real gamma);

    /**
     * @brief Diffusion discretization in 2D using central differences
     *
     * @param[in] A data to be discretized
     * @param[in] i x index
     * @param[in] j y index
     *
     */
    static Real diffusion(const Matrix<Real> &A, int i, int j);

    /**
     * @brief Convection in x direction using donor-cell scheme
     *
     * @param[in] U x-velocity field
     * @param[in] V y-velocity field
     * @param[in] i x index
     * @param[in] j y index
     * @param[out] result
     *
     */
    static Real convection_u(const Matrix<Real> &U, const Matrix<Real> &V, int i, int j);

    /**
     * @brief Convection in y direction using donor-cell scheme
     *
     * @param[in] U x-velocity field
     * @param[in] V y-velocity field
     * @param[in] i x index
     * @param[in] j y index
     * @param[out] result
     *
     */
    static Real convection_v(const Matrix<Real> &U, const Matrix<Real> &V, int i, int j);

    /**
     * @brief Convection of temperature in x direction using donor-cell scheme
     *
     * @param[in] U x-velocity field
     * @param[in] T Temperature
     * @param[in] i x index
     * @param[in] j y index
     * @param[out] result
     *
     */
    static Real convection_uT(const Matrix<Real> &U, const Matrix<Real> &T, int i, int j);

    /**
     * @brief Convection of temperature in y direction using donor-cell scheme
     *
     * @param[in] V y-velocity field
     * @param[in] T Temperature
     * @param[in] i x index
     * @param[in] j y index
     * @param[out] result
     *
     */
    static Real convection_vT(const Matrix<Real> &V, const Matrix<Real> &T, int i, int j);

    /**
     * @brief Convection of k and eps in x direction using donor-cell scheme
     *
     * @param[in] U x-velocity field
     * @param[in] KEPS k or epsilon
     * @param[in] i x index
     * @param[in] j y index
     * @param[out] result
     *
     */
    static Real convection_uKEPS(const Matrix<Real> &U, const Matrix<Real> &KEPS, int i, int j);

    /**
     * @brief Convection of k or epsilon in y direction using donor-cell scheme
     *
     * @param[in] V y-velocity field
     * @param[in] KEPS k or epsilon
     * @param[in] i x index
     * @param[in] j y index
     * @param[out] result
     *
     */
    static Real convection_vKEPS(const Matrix<Real> &V, const Matrix<Real> &KEPS, int i, int j);

    /**
     * @brief Laplacian term discretization for turbulent viscosity
     *
     * @param[in] P nu_t field
     * @param[in] nu molecular viscosity
     * @param[in] nu_i viscosity on vertical edges
     * @param[in] nu_j viscosity on horizontal edges
     * @param[in] i x index
     * @param[in] j y index
     * @param[in] coeff additional coefficient to multiply the result
     * @param[out] result
     *
     */
    static Real laplacian_nu(const Matrix<Real> &P, Real nu, const Matrix<Real> &nu_i, const Matrix<Real> &nu_j, int i,
                             int j, Real coeff = 1);

    /**
     * @brief Discretization of the mean strain rate term squared
     *
     * @param[in] U x-velocity field
     * @param[in] V y-velocity field
     * @param[in] S Shear stress field
     * @param[in] i x index
     * @param[in] j y index
     * @param[out] result
     *
     */
    static Real mean_strain_rate_squared(const Matrix<Real> &U, const Matrix<Real> &V, Matrix<Real> &S, int i, int j);

    /**
     * @brief Laplacian term discretization using central difference
     *
     * @param[in] P data to be discretized
     * @param[in] i x index
     * @param[in] j y index
     * @param[out] result
     *
     */
    static Real laplacian(const Matrix<Real> &P, int i, int j);

    /**
     * @brief Terms of laplacian needed for SOR, i.e. excluding unknown value at
     * (i,j)
     *
     * @param[in] P data to be discretized
     * @param[in] i x index
     * @param[in] j y index
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
     * @param[in] A data
     * @param[in] i x index
     * @param[in] j y index
     * @param[in] i_offset x offset
     * @param[in] j_offset y offset
     * @param[out] result
     */
    static Real diff(const Matrix<Real> &A, int i, int j, int i_offset, int j_offset);

    /// gamma value used for donor cell scheme
    static Real _gamma;

  private:
    /// cell width
    static Real _dx;
    /// cell height
    static Real _dy;
   
};
