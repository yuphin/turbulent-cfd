#include "Utilities.hpp"
#include <mpi.h>
#define LOG 1
class Communication {
  public:
    static void init_mpi(Params *params) {
        MPI_Init(NULL, NULL);
        MPI_Comm_size(MPI_COMM_WORLD, &params->world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &params->world_rank);
    }
    static void init_params(Params *params, int xmax, int jmax) {
        int local_size_x = xmax / params->world_size;
        int local_size_y = jmax / params->world_size;
        int iproc = params->iproc;
        int jproc = params->jproc;
        int wr = params->world_rank;
        Coord coords = {wr % params->iproc, wr / params->iproc};
        params->size_x = local_size_x;
        params->size_y = local_size_y;
        params->xmin = coords.x * local_size_x;
        params->xmax = params->xmin + local_size_x + 1;
        params->ymin = coords.y * local_size_y;
        params->ymax = params->ymin + local_size_y + 1;
        params->neighbor_left = coords.x == 0 ? -1 : wr - 1;
        params->neighbor_right = coords.x == iproc - 1 ? -1 : wr + 1;
        params->neighbor_bottom = coords.y == 0 ? -1 : wr - params->iproc;
        params->neighbor_top = coords.y == jproc - 1 ? -1 : wr + params->iproc;

#if LOG
        std::cout << "Rank: " << wr << "(" << coords.x << "," << coords.y << ")"
                  << "\n"
                  << "xmin: " << params->xmin << ", xmax: " << params->xmax << " ymin: " << params->ymin
                  << " ymax: " << params->ymax << "(Including ghost cells)"
                  << "\n"
                  << "Neighbors: \n"
                  << "Left: " << params->neighbor_left << " Right: " << params->neighbor_right
                  << " Top: " << params->neighbor_top << " Bottom: " << params->neighbor_bottom << "\n"
                  << "===================================" << std::endl;
#endif
    }
    static void finalize() { MPI_Finalize(); }
    static void communicate() {}
    static void reduce_min() {}
    static void reduce_sum() {}
};
