#pragma once

#include <array>
#include <vector>

#include "Utilities.hpp"

/**
 * @brief Main element of the grid that holds neighboring and type
 * information.
 *
 */
class Cell {
  public:
    Cell() = default;
    Cell(const Cell &other_cell) = default;
    /**
     * @brief Constructor for Cell object
     *
     * @param[in] x index of the cell
     * @param[in] y index of the cell
     * @param[in] type of the cell
     */
    Cell(int i, int j, cell_type type);

    /**
     * @brief Constructor for Cell object
     *
     * @param[in] x index of the cell
     * @param[in] y index of the cell
     * @param[in] type of the cell
     * @param[in] id of the cell, only for walls
     */
    Cell(int i, int j, cell_type type, int id);

    /**
     * @brief Neighbour getter for given border positon
     *
     * @param[in] border position
     * @param[out] pointer to the neighboring cell
     */
    const Cell *neighbour(border_position position) const;
    /**
     * @brief Neighbour setter for given border positon
     *
     * @param[in] border position
     */
    void set_neighbour(Cell *cell, border_position position);

    /**
     * @brief Access to the borders where fluid cells exist
     *
     * @param[out] vector of borders where fluid cells exist
     */
    const std::vector<border_position> &borders() const;

    /**
     * @brief Append new border to the cell
     *
     * Sets the corresponding border to true and adds the border
     * to the container
     *
     * @param[in] border position where fluid cell exists
     */
    void add_border(border_position border);
    /**
     * @brief Check whether the given position is a border to a
     * fluid cell or not
     *
     * @param[in] border position where fluid cell exists
     * @param[out] whether the given position is border or not
     */
    bool is_border(border_position position) const;

    /// Getter of x index
    int i() const;
    /// Getter of y index
    int j() const;
    /// Getter of cell type
    cell_type type() const;
    /// Getter of cell id
    int id() const;

    Real closest_dist;
    Real closest_wall_idx;

  private:
    /// x index
    int _i{0};
    /// y index
    int _j{0};
    /// Cell type
    cell_type _type{cell_type::DEFAULT};
    /// Cell id (only necessary for walls)
    int _id{0};

    /// Vector of bools that holds border conditions. TOP-BOTTOM-LEFT-RIGHT
    std::array<bool, 4> _border{false, false, false, false};
    /// Vector of border positions that holds existing borders
    std::vector<border_position> _borders;
    /// Pointers to neighbours. // TOP -  BOTTOM - LEFT - RIGHT - NORTHWEST -
    /// SOUTEAST
    std::array<Cell *, 6> _neighbours;
};
