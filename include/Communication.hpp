#pragma once

#include "Utilities.hpp"
#include <mpi.h>

#define LOG 1

class Communication {

  public:
    static void init_mpi(Params *params);
    static void init_params(Params *params, int imax, int jmax);
    
    static void communicate();
    static void reduce_min();
    static Real reduce_all(Real loc_value, MPI_Op mpi_operation);

    static void finalize();
};
