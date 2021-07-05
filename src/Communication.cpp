#include "Communication.hpp"

void Communication::init_mpi(Params *params) {
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &params->world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &params->world_rank);
}

void Communication::init_params(Params *params, int imax, int jmax) {
    int iproc = params->iproc;
    int jproc = params->jproc;    
    int wr = params->world_rank;
    Coord coords = {wr % params->iproc, wr / params->iproc};      
    int local_size_x = imax / iproc;
    int local_size_y = jmax / jproc;

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

    if (coords.x == iproc - 1) {
        params->size_x += imax % iproc;
        params->imax += imax % iproc;
    }
    if (coords.y == jproc - 1) {
        params->size_y += jmax % jproc;
        params->jmax += jmax % jproc;
    }   

#if FLUICHEN_LOG
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

void Communication::communicate(Params *params, Matrix<Real> &matrix) {
    int size_x = params->size_x;
    int size_y = params->size_y;

    MPI_Datatype column_type;
    int block_count_col = size_y + 2;
    int block_length_col = 1;
    int stride_col = size_x + 2;
    // TODO: Might need to use typedef later to adapt to different Real datatype.
    MPI_Type_vector(block_count_col, block_length_col, stride_col, MPI_DOUBLE, &column_type);
    MPI_Type_commit(&column_type);

    MPI_Datatype row_type;
    int block_count_row = size_x + 2;
    int block_length_row = 1;
    int stride_row = 1;
    // TODO: Might need to use typedef later to adapt to different Real datatype.
    MPI_Type_vector(block_count_row, block_length_row, stride_row, MPI_DOUBLE, &row_type);
    MPI_Type_commit(&row_type);

    // send left column
    MPI_Send(&matrix(1, 0), 1, column_type, params->neighbor_left, 0, MPI_COMM_WORLD);
    // receive right column
    MPI_Recv(&matrix(size_x + 1, 0), 1, column_type, params->neighbor_right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // send right column
    MPI_Send(&matrix(size_x, 0), 1, column_type, params->neighbor_right, 1, MPI_COMM_WORLD);
    // receive left column
    MPI_Recv(&matrix(0, 0), 1, column_type, params->neighbor_left, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // send top row
    MPI_Send(&matrix(0, size_y), 1, row_type, params->neighbor_top, 2, MPI_COMM_WORLD);
    // receive bottom row
    MPI_Recv(&matrix(0, 0), 1, row_type, params->neighbor_bottom, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // send bottom row
    MPI_Send(&matrix(0, 1), 1, row_type, params->neighbor_bottom, 3, MPI_COMM_WORLD);
    // receive top row
    MPI_Recv(&matrix(0, size_y + 1), 1, row_type, params->neighbor_top, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Type_free(&row_type);
    MPI_Type_free(&column_type);
}

Real Communication::reduce_all(Real loc_value, MPI_Op mpi_operation) {
    int rank;
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //Real global_value = 0.0;
    //// TODO: Might need to use typedef later to adapt to different Real datatype.
    //MPI_Allreduce(&loc_value, &global_value, 1, MPI_DOUBLE, mpi_operation, MPI_COMM_WORLD);
    return loc_value;
}

void Communication::finalize() { MPI_Finalize(); }