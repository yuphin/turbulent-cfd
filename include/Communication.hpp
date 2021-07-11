#pragma once

#include "Utilities.hpp"
#include "Datastructures.hpp"
#include <mpi.h>

#define FLUICHEN_LOG 0

class Communication {

  public:
    static void init_mpi(Params *params);
    static void init_params(Params *params, int imax, int jmax, const std::vector<Real> &dx_global, const std::vector<Real> &dy_global);
    
    static void communicate(Params *params, Matrix<Real> &matrix);
    static void reduce_min();
    static Real reduce_all(Real loc_value, MPI_Op mpi_operation);

    static void finalize();
};
