#pragma once

#include "Utilities.hpp"
#include "Datastructures.hpp"
#include <mpi.h>

#define FLUICHEN_LOG 0

/**
 * @brief Static Communication class that handles all MPI calls
 *
 */
class Communication {

  public:
    /**
     * @brief Initialize mpi
     * 
     * Set rank and worldsize
     * 
     * @param[in] params struct that holds mpi information
     */
    static void init_mpi(Params *params);

    /**
     * @brief Initialize the params struct with information needed for domain splitting
     *
     * @param[in] params struct that holds mpi information
     * @param[in] imax size in x of global domain
     * @param[in] jmax size in y of global domain
     */
    static void init_params(Params *params, int imax, int jmax);

    /**
     * @brief Communicate ghost layer values of matrix
     *
     * @param[in] params struct that holds mpi information
     * @param[in] matrix to communicate
     */
    static void communicate(Params *params, Matrix<Real> &matrix);

    /**
     * @brief Reduce values with specific operation
     *
     * @param[in] loc_value to be reduced
     * @param[in] mpi_operation to be used
     * @param[out] Real reduced value
     */
    static Real reduce_all(Real loc_value, MPI_Op mpi_operation);

    /**
     * @brief Finalize MPI
     */
    static void finalize();
};
