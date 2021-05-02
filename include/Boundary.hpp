#pragma once

#include "Cell.hpp"
#include "Fields.hpp"
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
    virtual void enforce_uv(Fields &field) = 0;
    virtual void enforce_fg(Fields &field) = 0;
    virtual void enforce_p(Fields &field) = 0;
    virtual ~Boundary() = default;
};

/**
 * @brief Fixed wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, which is zero, Neumann for pressure
 */
class FixedWallBoundary : public Boundary {
  public:
    FixedWallBoundary(std::vector<Cell *> cells);
    FixedWallBoundary(std::vector<Cell *> cells, std::unordered_map<int, double> wall_temperature);
    virtual ~FixedWallBoundary() = default;
    void enforce_uv(Fields &field) override;
    void enforce_fg(Fields &field) override;
    void enforce_p(Fields &field) override;

  private:
    std::vector<Cell *> _cells;
    std::unordered_map<int, double> _wall_temperature;
};

/**
 * @brief Moving wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities for the given velocity parallel to the fluid,
 * Neumann for pressure
 */
class MovingWallBoundary : public Boundary {
  public:
    MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity);
    MovingWallBoundary(std::vector<Cell *> cells, std::unordered_map<int, double> wall_velocity,
                       std::unordered_map<int, double> wall_temperature);
    virtual ~MovingWallBoundary() = default;
    void enforce_uv(Fields &field) override;
    void enforce_fg(Fields &field) override;
    void enforce_p(Fields &field) override;    

  private:
    std::vector<Cell *> _cells;
    std::unordered_map<int, double> _wall_velocity;
    std::unordered_map<int, double> _wall_temperature;
};
