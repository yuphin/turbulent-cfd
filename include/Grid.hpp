#pragma once

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Cell.hpp"
#include "Datastructures.hpp"
#include "Domain.hpp"
#include "Utilities.hpp"

/**
 * @brief Data structure holds cells and related sub-containers
 *
 */
class Grid {
  public:
    Grid() = default;

    /**
     * @brief Constructor for the Grid
     *
     * @param[in] domain object specifying the domain
     * @param[in] geometry_data layout of the problem case
     *
     */
    Grid(Domain &domain, std::vector<std::vector<int>> &geometry_data);

    /// index based cell access
    Cell cell(int i, int j) const;

    /// access number of cells in x direction
    int imax() const;
    /// access number of cells in y direction
    int jmax() const;
    /// access number of cells in x direction including ghost cells
    int imaxb() const;
    /// access number of cells in x direction including ghost cells
    int jmaxb() const;

    /// access number of cells in x direction excluding ghost cells
    const Domain &domain() const;

    /// access cell size in x-direction
    Real dx() const;
    /// access cell size in y-direction
    Real dy() const;

    /**
     * @brief Access inflow cells
     *
     * @param[out] vector of fluid cells
     */
    const std::vector<Cell *> &fluid_cells() const;

    /**
     * @brief Access inlet cells
     *
     * @param[out] vector of inlet cells
     */
    std::vector<Cell *> &inlet_cells();

    /**
     * @brief Access outlet cells
     *
     * @param[out] vector of outlet cells
     */
    std::vector<Cell *> &outlet_cells();

    /**
     * @brief Access noslip wall cells
     *
     * @param[out] vector of moving wall cells
     */
    std::vector<Cell *> &noslip_wall_cells();

    /**
     * @brief Access freeslip cells
     *
     * @param[out] vector of fixed wall cells
     */
    std::vector<Cell *> &freeslip_wall_cells();

    /// domain object 
    Domain _domain;

    /// All cells of the grid
    Matrix<Cell> _cells;

    /// Calculate cell wall distances
    void preprocess_geometry();

  private:
    /// Build cell data structures with given geometrical data
    void assign_cell_types(std::vector<std::vector<int>> &geometry_data);

    /// Cell pointer to all fluid cells
    std::vector<Cell *> _fluid_cells;
    /// Cell pointer to all inlet cells
    std::vector<Cell *> _inlet_cells;
    /// Cell pointer to all outlet cells
    std::vector<Cell *> _outlet_cells;
    /// Cell pointer to all freeslip wall cells
    std::vector<Cell *> _freeslip_wall_cells;
    /// Cell pointer to all noslip wall cells
    std::vector<Cell *> _noslip_wall_cells;

    /// Cells width
    Real _dx;
    /// Cells height
    Real _dy;
};
