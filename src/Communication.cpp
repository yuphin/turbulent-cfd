#include "Communication.hpp"


void Communication::init_mpi(Params *params) {
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &params->world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &params->world_rank);
}

void Communication::init_params(Params *params, int imax, int jmax) {
        int iproc = params->iproc;
        int jproc = params->jproc;
        int local_size_x = imax / iproc;
        int local_size_y = jmax / jproc;
        int wr = params->world_rank;
        Coord coords = {wr % params->iproc, wr / params->iproc};
        params->size_x = local_size_x;
        params->size_y = local_size_y;
        params->imin = coords.x * local_size_x;
        params->imax = params->imin + local_size_x + 1;
        params->jmin = coords.y * local_size_y;
        params->jmax = params->jmin + local_size_y + 1;
        params->neighbor_left = coords.x == 0 ? MPI_PROC_NULL : wr - 1;
        params->neighbor_right = coords.x == iproc - 1 ? MPI_PROC_NULL : wr + 1;
        params->neighbor_bottom = coords.y == 0 ? MPI_PROC_NULL : wr - params->iproc;
        params->neighbor_top = coords.y == jproc - 1 ? MPI_PROC_NULL : wr + params->iproc;

#if LOG
        std::cout << "Rank: " << wr << "(" << coords.x << "," << coords.y << ")"
                  << "\n"
                  << "xmin: " << params->imin << ", xmax: " << params->imax << " ymin: " << params->jmin
                  << " ymax: " << params->jmax << "(Including ghost cells)"
                  << "\n"
                  << "Neighbors: \n"
                  << "Left: " << params->neighbor_left << " Right: " << params->neighbor_right
                  << " Top: " << params->neighbor_top << " Bottom: " << params->neighbor_bottom << "\n"
                  << "===================================" << std::endl;
#endif
}


Real Communication::reduce_all(Real loc_value, MPI_Op mpi_operation) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Real global_value = 0.0;
    MPI_Allreduce(&loc_value, &global_value, 1, MPI_DOUBLE, mpi_operation, MPI_COMM_WORLD);
    // std::cout << "\nRank: " << rank << "   " << global_value << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    return global_value;
}

void Communication::finalize() { MPI_Finalize(); }