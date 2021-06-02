#pragma once

#include "Utilities.hpp"
#define LOG 1

class Communication {

  public:
    static void init_mpi(Params *params);
    static void init_params(Params *params, int imax, int jmax);
    
    static void communicate();
    static void reduce_min();
    static void reduce_sum();

    static void finalize();
};
