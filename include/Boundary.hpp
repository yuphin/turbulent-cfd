#pragma once

#include "Cell.hpp"
#include "Fields.hpp"
#include <unordered_map>
#include <vector>
#include <cassert>

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
    virtual void enforce_uv(Fields &field) = 0;
    virtual void enforce_fg(Fields &field);
    virtual void enforce_p(Fields &field);
    virtual void enforce_t(Fields &field);
    virtual ~Boundary() = default;

  protected:
    virtual void enforce_p1(Fields &field, Cell* cell);
    virtual void enforce_p2(Fields &field, Cell* cell);
    std::vector<Cell *> *_cells;
    std::unordered_map<int, double> _wall_temperature;
  private:
};

/**
 * @brief Fixed wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, which is zero, Neumann for pressure
 */
class OutletBoundary : public Boundary {
  public:
    OutletBoundary(std::vector<Cell *> *cells);
    virtual ~OutletBoundary() = default;
    void enforce_uv(Fields &field) override;
    void enforce_t(Fields &field) override { assert(false); }
};

/**
 * @brief Fixed wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, which is zero, Neumann for pressure
 */
class InletBoundary : public Boundary {
  public:
    InletBoundary(std::vector<Cell *> *cells);
    InletBoundary(std::vector<Cell *> *cells, std::unordered_map<int, double> inlet_U,
                      std::unordered_map<int, double> inlet_V,
                      std::unordered_map<int, double> inlet_T);
    virtual ~InletBoundary() = default;
    void enforce_uv(Fields &field) override;
    void enforce_t(Fields &field) override { assert(false); }

  private:
    std::unordered_map<int, double> _inlet_U;
    std::unordered_map<int, double> _inlet_V;
    std::unordered_map<int, double> _inlet_T;    
};


/**
 * @brief Fixed wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, which is zero, Neumann for pressure
 */
class NoSlipWallBoundary : public Boundary {
  public:
    NoSlipWallBoundary(std::vector<Cell *> *cells);
    NoSlipWallBoundary(std::vector<Cell *> *cells, std::unordered_map<int, double> wall_velocity,
                      std::unordered_map<int, double> wall_temperature);
    virtual ~NoSlipWallBoundary() = default;
    void enforce_uv(Fields &field) override;
   
  private:
    void enforce_uv1(Fields &field, Cell* cell);
    void enforce_uv2(Fields &field, Cell* cell);
    std::unordered_map<int, double> _wall_velocity;
};

/**
 * @brief Moving wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities for the given velocity parallel to the fluid,
 * Neumann for pressure
 */
class FreeSlipWallBoundary : public Boundary {
  public:
    FreeSlipWallBoundary(std::vector<Cell *> *cells);
    FreeSlipWallBoundary(std::vector<Cell *> *cells, std::unordered_map<int, double> wall_velocity,
                       std::unordered_map<int, double> wall_temperature);
    virtual ~FreeSlipWallBoundary() = default;
    void enforce_uv(Fields &field) override;

  private:
    std::unordered_map<int, double> _wall_velocity;
};
