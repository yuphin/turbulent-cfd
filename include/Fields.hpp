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
    Fields(double _nu, double _dt, double _tau, int imax, int jmax, 
           double UI, double VI, double PI, double TI, double _alpha, double _beta,
           double _gx, double _gy);

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
    double calculate_dt(Grid &grid, bool calc_temp);

    /// x-velocity index based access and modify
    double &u(int i, int j);

    /// y-velocity index based access and modify
    double &v(int i, int j);

    /// pressre index based access and modify
    double &p(int i, int j);

    /// temperature index based acccess and modify
    double &t(int i, int j);

    /// RHS index based access and modify
    double &rs(int i, int j);

    /// x-momentum flux index based access and modify
    double &f(int i, int j);

    /// y-momentum flux index based access and modify
    double &g(int i, int j);

    /// get timestep size
    double dt() const;

    /// initial pressure
    double PI;

    /// initial temperature
    double TI;

    /// pressure matrix access and modify
    Matrix<double> &p_matrix();

    /// velocity u matrix access and modify
    Matrix<double> &u_matrix();

    /// velocity v matrix access and modify
    Matrix<double> &v_matrix();

    /// temperature t matrix access and modify
    Matrix<double> &t_matrix();

    /// x-momentum flux matrix
    Matrix<double> &f_matrix();

    /// y-momentum flux matrix
    Matrix<double> &g_matrix();

  private:
    /// x-velocity matrix
    Matrix<double> _U;
    /// y-velocity matrix
    Matrix<double> _V;
    /// pressure matrix
    Matrix<double> _P;
    /// temperature matrix
    Matrix<double> _T;
    /// x-momentum flux matrix
    Matrix<double> _F;
    /// y-momentum flux matrix
    Matrix<double> _G;
    /// right hand side matrix
    Matrix<double> _RS;

    /// kinematic viscosity
    double _nu;
    /// gravitional accelearation in x direction
    double _gx{0.0};
    /// gravitional accelearation in y direction
    double _gy{0.0};
    /// timestep size
    double _dt;
    /// adaptive timestep coefficient
    double _tau;
    /// thermal expansion coefficient
    double _beta;
    /// thermal diffusivity
    double _alpha;
};
