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
     * @param[in] geometry file name
     * @param[in] number of cells in x direction
     * @param[in] number of cells in y direction
     * @param[in] cell size in x direction
     * @param[in] cell size in y direction
     *
     */
    Grid(std::string geom_name, Domain &domain, 
         std::vector<std::vector<int>> &geometry_data);

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
     * @brief Access moving wall cells
     *
     * @param[out] vector of moving wall cells
     */
    std::vector<Cell *> &noslip_wall_cells();

    /**
     * @brief Access fixed wall cells
     *
     * @param[out] vector of fixed wall cells
     */
    std::vector<Cell *> &freeslip_wall_cells();

    Domain _domain;
    Matrix<Cell> _cells;

     void preprocess_geometry() {
        auto at = [this](int i, int j) { return &_cells._container[j * _domain.imax + i]; };
        auto get_idx = [this](int i, int j) { return j * _domain.imax + i; };
        for (auto cell : _fluid_cells) {
            int i = cell->i();
            int j = cell->j();
            cell_type ct;
            Real dist_x_p = 0;
            Real dist_y_p = 0;
            Real dist_x_n = 0;
            Real dist_y_n = 0;
            int idx_i_p = i;
            do {
                idx_i_p++;
                ct = at(idx_i_p, j)->type();
                // dist_x_p += at(idx_i_p, j)->dx;
                dist_x_p += dx();
            } while (ct == cell_type::FLUID);
            if (ct != cell_type::FREESLIP_WALL && ct != cell_type::NOSLIP_WALL) {
                dist_x_p = REAL_MAX;
            }
            int idx_i_n = i;
            do {
                idx_i_n--;
                ct = at(idx_i_n, j)->type();
                // dist_x_n += at(idx_i_n, j)->dx;
                dist_x_n += dx();
            } while (ct == cell_type::FLUID);
            if (ct != cell_type::FREESLIP_WALL && ct != cell_type::NOSLIP_WALL) {
                dist_x_n = REAL_MAX;
            }
            int idx_j_p = j;
            do {
                idx_j_p++;
                ct = at(i, idx_j_p)->type();
                // dist_y_p += at(i, idx_j_p)->dy;
                dist_y_p += dy();
            } while (ct == cell_type::FLUID);
            if (ct != cell_type::FREESLIP_WALL && ct != cell_type::NOSLIP_WALL) {
                dist_y_p = REAL_MAX;
            }
            int idx_j_n = j;
            do {
                idx_j_n--;
                ct = at(i, idx_j_n)->type();
                // dist_y_n += at(i, idx_j_n)->dy;
                dist_y_n += dy();
            } while (ct == cell_type::FLUID);
            if (ct != cell_type::FREESLIP_WALL && ct != cell_type::NOSLIP_WALL) {
                dist_y_n = REAL_MAX;
            }
            Cell *cell = at(i, j);
            std::array<int, 4> idxs = {get_idx(idx_i_p, j), get_idx(idx_i_n, j), get_idx(i, idx_j_p),
                                       get_idx(i, idx_j_n)};
            std::array<Real, 4> dists = {dist_x_p, dist_x_n, dist_y_p, dist_y_n};
            Real *min_xy = std::min_element(dists.data(), dists.data() + 4);
            int min_idx = min_xy - dists.data();
            cell->closest_dist = *min_xy;
            cell->closest_wall_idx = idxs[min_idx];
        }
    }

  private:
    /// Build cell data structures with given geometrical data
    void assign_cell_types(std::vector<std::vector<int>> &geometry_data);

    std::vector<Cell *> _fluid_cells;
    std::vector<Cell *> _inlet_cells;
    std::vector<Cell *> _outlet_cells;
    std::vector<Cell *> _freeslip_wall_cells;
    std::vector<Cell *> _noslip_wall_cells;


    Real _dx;
    Real _dy;
};
