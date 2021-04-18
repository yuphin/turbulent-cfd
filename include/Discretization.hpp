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
    Discretization(double dx, double dy, double gamma);

    /**
     * @brief Diffusion discretization in 2D using central differences
     *
     * @param[in] data to be discretized
     * @param[in] x index
     * @param[in] y index
     *
     */
    static double diffusion(const Matrix<double> &A, int i, int j);

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
    static double convection_u(const Matrix<double> &U, const Matrix<double> &V, int i, int j);

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
    static double convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j);

    /**
     * @brief Laplacian term discretization using central difference
     *
     * @param[in] data to be discretized
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static double laplacian(const Matrix<double> &P, int i, int j);

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
    static double sor_helper(const Matrix<double> &P, int i, int j);

    /**
     * @brief Linear interpolation
     *
     * @param[in] data to be interpolated
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static double interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset);

  private:
    static double _dx;
    static double _dy;
    static double _gamma;
};
