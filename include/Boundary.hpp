#pragma once
#include "Cell.hpp"
#include "Fields.hpp"
#include "Utilities.hpp"
#include <cassert>
#include <unordered_map>
#include <vector>

/**
 * @brief Abstact of boundary conditions.
 *
 * This class patches the physical values to the given field.
 * Default Boundary Conditions are specified here.
 */
class Boundary {
  public:
    /**
     * @brief Constructor for a Boundary object
     *
     * @param[in] Pointer to a vector of Cell pointers
     */
    Boundary(std::vector<Cell *> *cells) : _cells(cells){};

    /**
     * @brief Default method to handle u and v boundary condition
     *
     * u and v are not handled by default
     * 
     * @param[in] Field to be applied
     */
    virtual void enforce_uv(Fields &field){};

    /**
     * @brief Default method to handle f and g boundary condition
     *
     * Set f and g values to corresponding u and v values
     * 
     * @param[in] Field to be applied
     */
    virtual void enforce_fg(Fields &field);

    /**
     * @brief Default method to handle pressure boundary condition
     *
     * Set zero von Neumann condition
     * 
     * @param[in] Field to be applied
     */
    virtual void enforce_p(Fields &field);

    /**
     * @brief Default method to handle temperature boundary condition
     *
     * Set Dirichlet condition to _wall_temperature
     * or adiabatic condition if _wall_temperature[cell->id] == -1
     * 
     * @param[in] Field to be applied
     */
    virtual void enforce_t(Fields &field);

    /**
     * @brief Default method to handle turbulence boundary condition
     *
     * @param[in] Field to be applied
     * @param[in] Turbulence model that is being used
     */
    virtual void enforce_nu_t(Fields &field, int turb_model);

    /// Getter of boundary type 
    uint32_t get_type() { return _type; }

    /// Vector pointer holding the cells of the boundary
    std::vector<Cell *> *_cells;
    /// Map holding wall temperatures
    std::unordered_map<int, Real> _wall_temperature;
    /// value holding the type of boundary
    uint32_t _type;
    virtual ~Boundary() = default;

  protected:
    /**
     * @brief Apply von Neumann pressure boundary condition on
     *        cell with one neighboring fluid cell
     *
     * @param[in] Field to be applied
     * @param[i] Cell pointer to be applied
     */
    virtual void enforce_p_main(Fields &field, Cell *cell);

    /**
     * @brief Apply von Neumann pressure boundary condition on
     *        cell with two neighboring fluid cells
     *
     * @param[in] Field to be applied
     * @param[i] Cell pointer to be applied
     */
    virtual void enforce_p_diagonal(Fields &field, Cell *cell);

    /**
     * @brief Apply Dirichlet temperature boundary condition on
     *        cell with one neighboring fluid cell
     *
     * @param[in] Field to be applied
     * @param[i] Cell pointer to be applied
     */
    void enforce_t_drichlet_main(Fields &field, Cell *cell);

    /**
     * @brief Apply Dirichlet temperature boundary condition on
     *        cell with two neighboring fluid cells
     *
     * @param[in] Field to be applied
     * @param[i] Cell pointer to be applied
     */
    void enforce_t_drichlet_diag(Fields &field, Cell *cell);

    /**
     * @brief Apply adiabatic temperature boundary condition on
     *        cell with one neighboring fluid cell
     *
     * @param[in] Field to be applied
     * @param[i] Cell pointer to be applied
     */
    void enforce_t_adiabatic_main(Fields &field, Cell *cell);

    /**
     * @brief Apply adiabatic temperature boundary condition on
     *        cell with one neighboring fluid cell
     *
     * @param[in] Field to be applied
     * @param[i] Cell pointer to be applied
     */
    void enforce_t_adiabatic_diag(Fields &field, Cell *cell);
};

/**
 * @brief Outlet boundary condition where fluid exits the domain.
 * Neumann for velocities, Neumann for pressure, Neumann for temperature
 */
class OutletBoundary : public Boundary {
  public:
    /**
     * @brief Constructor for an Outlet boundary object
     *
     * @param[in] Pointer to a vector of Cell pointers
     */
    OutletBoundary(std::vector<Cell *> *cells);

    /**
     * @brief Apply Dirichlet pressure condition setting value to initial pressure
     *
     * @param[in] Field to be applied
     */
    void enforce_p(Fields &field) override;

    /**
     * @brief Apply adiabatic temperature condition
     *
     * @param[in] Field to be applied
     */
    void enforce_t(Fields &field) override;
    virtual ~OutletBoundary() = default;    
};

/**
 * @brief Inlet boundary condition where fluid enters the domain.
 * Dirichlet for velocities, Neumann for pressure, Dirichlet for temperature
 */
class InletBoundary : public Boundary {
  public:
    InletBoundary(std::vector<Cell *> *cells);

    /**
     * @brief Constructor for InletBoundary object
     *
     * @param[in] cells vector pointer to Cell pointers
     * @param[in] inlet_U map of inlet u velocities
     * @param[in] inlet_V map of inlet v velocities
     * @param[in] inlet_T map of inlet temperatures
     * @param[in] inlet_K map of inlet k values
     * @param[in] inlet_EPS map of inlet epsilon values
     * @param[in] DP inlet pressure value
     */
    InletBoundary(std::vector<Cell *> *cells, std::unordered_map<int, Real> inlet_U,
                  std::unordered_map<int, Real> inlet_V, std::unordered_map<int, Real> inlet_T,
                  std::unordered_map<int, Real> inlet_K, std::unordered_map<int, Real> inlet_EPS, Real DP);
    virtual ~InletBoundary() = default;

    /**
     * @brief Apply fixed input velocities
     *
     * @param[in] Field to be applied
     */
    void enforce_uv(Fields &field) override;

    /**
     * @brief Apply Dirichlet pressure boundary condition if _inlet_DP is set
     *
     * @param[in] Field to be applied
     */
    void enforce_p(Fields &field) override;

    /**
     * @brief Apply Dirichlet temperature boundary condition
     *
     * @param[in] Field to be applied
     */
    void enforce_t(Fields &field) override;

    /**
     * @brief Apply Dirichlet turbulence boundary conditions
     *
     * @param[in] Field to be applied
     * @param[in] turb_model to be used
     */
    void enforce_nu_t(Fields &field, int turb_model) override;

    /// Map holding u velocities of all inlets
    std::unordered_map<int, Real> _inlet_U;
    /// Map holding v velocities of all inlets
    std::unordered_map<int, Real> _inlet_V;
    /// Map holding temperatures of all inlets
    std::unordered_map<int, Real> _inlet_T;
    /// Map holding k values of all inlets
    std::unordered_map<int, Real> _inlet_K;
    /// Map holding epsilon values of all inlets
    std::unordered_map<int, Real> _inlet_EPS;

  private:
    /// Pressure value if using pressure inlet
    Real _inlet_DP = REAL_MAX;
};

/**
 * @brief Fixed wall boundary condition for the boundaries of the domain.
 * Dirichlet for velocities, which is zero or the wall velocity, Neumann for pressure
 */
class NoSlipWallBoundary : public Boundary {
  public:
    NoSlipWallBoundary(std::vector<Cell *> *cells);

    /**
     * @brief Constructor for NoSlipBoundary
     *
     * @param[in] cells vector pointer to Cell pointers
     * @param[in] wall_velocity Map with tangential wall velocities
     * @param[in] wall_temperature Map with wall temperatures
     */
    NoSlipWallBoundary(std::vector<Cell *> *cells, std::unordered_map<int, Real> wall_velocity,
                       std::unordered_map<int, Real> wall_temperature);

    /**
     * @brief Apply no-slip boundary condition
     *
     * Set normal velocities 0 and tangential to wall_velocity
     * 
     * @param[in] Fiels to be applied
     */    
    void enforce_uv(Fields &field) override;

    /// tangential wall velocities
    std::unordered_map<int, Real> _wall_velocity;
    virtual ~NoSlipWallBoundary() = default;

  private:
    /**
     * @brief Apply no-slip boundary condition for cell with
     *        one neigboring fluid cell
     *
     * @param[in] Fiels to be applied
     * @param[in] cell to be applid
     */  
    void enforce_uv_main(Fields &field, Cell *cell);

    /**
     * @brief Apply no-slip boundary condition for cell with
     *        two neigboring fluid cells
     *
     * @param[in] Fiels to be applied
     * @param[in] cell to be applied
     */  
    void enforce_uv_diagonal(Fields &field, Cell *cell);
};

/**
 * @brief Freeslip wall boundary condition for the boundaries of the domain.
 * Dirichlet for velocities for the given velocity parallel to the fluid,
 * Neumann for pressure
 */
class FreeSlipWallBoundary : public Boundary {
  public:
    FreeSlipWallBoundary(std::vector<Cell *> *cells);

    /**
     * @brief Constructor for FreeSlipWallBoundary
     *
     * @param[in] cells vector pointer to Cell pointers
     * @param[in] wall_velocity Map with tangential wall velocities
     * @param[in] wall_temperature Map with wall temperatures
     */
    FreeSlipWallBoundary(std::vector<Cell *> *cells, std::unordered_map<int, Real> wall_velocity,
                         std::unordered_map<int, Real> wall_temperature);
    virtual ~FreeSlipWallBoundary() = default;

    /**
     * @brief Apply zero Dirichlet boundary condition only in
     *        normal direction
     *
     * @param[in] Field to be applied
     */
    void enforce_uv(Fields &field) override;

    /// Map with tangential wall velocities
    std::unordered_map<int, Real> _wall_velocity;

};
