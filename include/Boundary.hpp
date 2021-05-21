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
 */
class Boundary {
  public:
    /**
     * @brief Main method to patch the boundary conditons to given field and
     * grid
     *
     * @param[in] Field to be applied
     */
    Boundary(std::vector<Cell *> *cells) : _cells(cells){};
    virtual void enforce_uv(Fields &field){}; // Don't handle by default
    virtual void enforce_fg(Fields &field);
    virtual void enforce_p(Fields &field);
    virtual void enforce_t(Fields &field);
    virtual ~Boundary() = default;
    std::vector<Cell *> *_cells;
    std::unordered_map<int, Real> _wall_temperature;

  protected:
    virtual void enforce_p_main(Fields &field, Cell *cell);
    virtual void enforce_p_diagonal(Fields &field, Cell *cell);

  private:
    void enforce_t_drichlet_main(Fields &field, Cell *cell);
    void enforce_t_drichlet_diag(Fields &field, Cell *cell);
    void enforce_t_outflow_main(Fields &field, Cell *cell);
    void enforce_t_outflow_diag(Fields &field, Cell *cell);
};

/**
 * @brief Fixed wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, which is zero, Neumann for pressure
 */
class OutletBoundary : public Boundary {
  public:
    OutletBoundary(std::vector<Cell *> *cells);
    virtual ~OutletBoundary() = default;
    void enforce_t(Fields &field) override {}
    void enforce_p(Fields &field) override;
};

/**
 * @brief Fixed wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, which is zero, Neumann for pressure
 */
class InletBoundary : public Boundary {
  public:
    InletBoundary(std::vector<Cell *> *cells);
    InletBoundary(std::vector<Cell *> *cells, std::unordered_map<int, Real> inlet_U,
                  std::unordered_map<int, Real> inlet_V, std::unordered_map<int, Real> inlet_T, Real DP);
    virtual ~InletBoundary() = default;
    void enforce_uv(Fields &field) override;
    void enforce_t(Fields &field) override {}
    void enforce_p(Fields &field) override;

  private:
    std::unordered_map<int, Real> _inlet_U;
    std::unordered_map<int, Real> _inlet_V;
    std::unordered_map<int, Real> _inlet_T;
    Real _inlet_DP = REAL_MAX;
};

/**
 * @brief Fixed wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, which is zero, Neumann for pressure
 */
class NoSlipWallBoundary : public Boundary {
  public:
    NoSlipWallBoundary(std::vector<Cell *> *cells);
    NoSlipWallBoundary(std::vector<Cell *> *cells, std::unordered_map<int, Real> wall_velocity,
                       std::unordered_map<int, Real> wall_temperature);
    virtual ~NoSlipWallBoundary() = default;
    void enforce_uv(Fields &field) override;

  private:
    void enforce_uv_main(Fields &field, Cell *cell);
    void enforce_uv_diagonal(Fields &field, Cell *cell);
    std::unordered_map<int, Real> _wall_velocity;
};

/**
 * @brief Moving wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities for the given velocity parallel to the fluid,
 * Neumann for pressure
 */
class FreeSlipWallBoundary : public Boundary {
  public:
    FreeSlipWallBoundary(std::vector<Cell *> *cells);
    FreeSlipWallBoundary(std::vector<Cell *> *cells, std::unordered_map<int, Real> wall_velocity,
                         std::unordered_map<int, Real> wall_temperature);
    virtual ~FreeSlipWallBoundary() = default;
    void enforce_uv(Fields &field) override;

  private:
    std::unordered_map<int, Real> _wall_velocity;
};
