#include "PressureSolver.hpp"
#include "Communication.hpp"

#include <cmath>
#include <iostream>
#include <mpi.h>

SOR::SOR(Real omega) : _omega(omega) {}

Real SOR::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries, Params &params) {

    Real dx = grid.dx();
    Real dy = grid.dy();

    Real coeff = _omega / (2.0f * (1.0 / (dx * dx) + 1.0 / (dy * dy))); // = _omega * h^2 / 4.0, if dx == dy == h

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        field.p(i, j) = (1.0 - _omega) * field.p(i, j) +
                        coeff * (Discretization::sor_helper(field.p_matrix(), i, j) - field.rs(i, j));
    }

    Real res = 0.0;
    Real rloc = 0.0;

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        Real val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
        rloc += (val * val);
    }
    // Exchange pressure
    Communication::communicate(&params, field.p_matrix());
    // Compute global residual
    int global_cells;
    res = Communication::reduce_all(rloc, MPI_SUM);
    global_cells = Communication::reduce_all(grid.fluid_cells().size(), MPI_SUM);
    {
        res = res / global_cells;
        res = std::sqrt(res);
    }
    
    return res;
}
