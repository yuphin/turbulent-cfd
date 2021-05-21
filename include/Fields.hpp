#pragma once

#include "Datastructures.hpp"
#include "Discretization.hpp"
#include "Grid.hpp"

/**
 * @brief Class of container and modifier for the physical fields
 *
 */
class Fields {
  public:
    Fields() = default;

    /**
     * @brief Constructor for the fields
     *
     * @param[in] kinematic viscosity
     * @param[in] initial timestep size
     * @param[in] adaptive timestep coefficient
     * @param[in] number of cells in x direction
     * @param[in] number of cells in y direction
     * @param[in] initial x-velocity
     * @param[in] initial y-velocity
     * @param[in] initial pressure
     *
     */
    Fields(Real _nu, Real _dt, Real _tau, int imax, int jmax, 
           Real UI, Real VI, Real PI, Real TI, Real _alpha, Real _beta,
           Real _gx, Real _gy);

    /**
     * @brief Calculates the convective and diffusive fluxes in x and y
     * direction based on explicit discretization of the momentum equations
     *
     * @param[in] grid in which the fluxes are calculated
     * @param[in] whether to include temperature related terms
     *
     */
    void calculate_fluxes(Grid &grid, bool calc_temp);

    /**
     * @brief Right hand side calculations using the fluxes for the pressure
     * Poisson equation
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_rs(Grid &grid);

    /**
     * @brief Velocity calculation using pressure values
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_velocities(Grid &grid);

    
    /**
     * @brief Velocity calculation using pressure values
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_temperatures(Grid &grid);

    /**
     * @brief Adaptive step size calculation using x-velocity condition,
     * y-velocity condition and CFL condition
     *
     * @param[in] grid in which the calculations are done
     * @param[in] whether to include temperatures
     */
    Real calculate_dt(Grid &grid, bool calc_temp);

    /// x-velocity index based access and modify
    Real &u(int i, int j);

    /// y-velocity index based access and modify
    Real &v(int i, int j);

    /// pressre index based access and modify
    Real &p(int i, int j);

    /// temperature index based acccess and modify
    Real &t(int i, int j);

    /// RHS index based access and modify
    Real &rs(int i, int j);

    /// x-momentum flux index based access and modify
    Real &f(int i, int j);

    /// y-momentum flux index based access and modify
    Real &g(int i, int j);

    /// get timestep size
    Real dt() const;

    /// initial pressure
    Real PI;

    /// initial temperature
    Real TI;

    /// pressure matrix access and modify
    Matrix<Real> &p_matrix();

    /// velocity u matrix access and modify
    Matrix<Real> &u_matrix();

    /// velocity v matrix access and modify
    Matrix<Real> &v_matrix();

    /// temperature t matrix access and modify
    Matrix<Real> &t_matrix();

    /// x-momentum flux matrix
    Matrix<Real> &f_matrix();

    /// y-momentum flux matrix
    Matrix<Real> &g_matrix();

  private:
    /// x-velocity matrix
    Matrix<Real> _U;
    /// y-velocity matrix
    Matrix<Real> _V;
    /// pressure matrix
    Matrix<Real> _P;
    /// temperature matrix
    Matrix<Real> _T;
    /// x-momentum flux matrix
    Matrix<Real> _F;
    /// y-momentum flux matrix
    Matrix<Real> _G;
    /// right hand side matrix
    Matrix<Real> _RS;

    /// kinematic viscosity
    Real _nu;
    /// gravitional accelearation in x direction
    Real _gx{0.0};
    /// gravitional accelearation in y direction
    Real _gy{0.0};
    /// timestep size
    Real _dt;
    /// adaptive timestep coefficient
    Real _tau;
    /// thermal expansion coefficient
    Real _beta;
    /// thermal diffusivity
    Real _alpha;
};
