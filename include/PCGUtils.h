#pragma once
#include <lib/pcgsolver.h>
#include "Grid.hpp"
#include "Boundary.hpp"

DiagonalSparseMatrix<Real> create_diagonal_matrix(const SparseMatrix<Real> &A, int dim_x, int dim_y,
                                                  const std::vector<int> offsets);

DiagonalSparseMatrix<Real> create_preconditioner_spai(const SparseMatrix<Real> &A, Grid &_grid, int precond_type = 0);
void build_pcg_matrix(Fields &_field, Grid &_grid, const std::vector<std::unique_ptr<Boundary>> &_boundaries,
                      SparseMatrix<Real> &A, SparseMatrix<Real> &U, SparseMatrix<Real> &V, SparseMatrix<Real> &T,
                      std::vector<Real> &U_RHS, std::vector<Real> &V_RHS, std::vector<Real> &T_RHS,
                      FixedSparseMatrix<Real> &U_fixed, FixedSparseMatrix<Real> &V_fixed,
                      FixedSparseMatrix<Real> &T_fixed);
