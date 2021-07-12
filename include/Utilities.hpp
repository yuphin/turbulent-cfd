#pragma once

#include <cfloat>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <mpi.h>

#define USE_FLOATS 0

#if USE_FLOATS
const float REAL_MAX = FLT_MAX;
typedef float Real;
#else
const double REAL_MAX = DBL_MAX;
typedef double Real;
#endif

// If no geometry file is provided in the input file, lid driven cavity case
// will run by default. In the Grid.cpp, geometry will be created following
// PGM convention, which is:
// 0: fluid, 10: fixed wall, 11: moving wall
namespace LidDrivenCavity {
const int moving_wall_id = 11;
const int fixed_wall_id = 10;
const Real wall_velocity = 1.0;
} // namespace LidDrivenCavity


//// ENUMS and STRUCTS ////

/// Solver type enum
enum class SolverType { SOR, PCG };

/// Border position enumm
enum class border_position {
    TOP,
    BOTTOM,
    LEFT,
    RIGHT,
};

/// Boundary type enum
enum BoundaryType {
    BND_NOSLIP_IDX,
    BND_FREESLIP_IDX,
    BND_OUTLET_IDX,
    BND_INLET_IDX
};

/// Descriptive border definitions
namespace border {
    const int TOP = 0;
    const int BOTTOM = 1;
    const int LEFT = 2;
    const int RIGHT = 3;
} // namespace border

/// Cell type enum
enum class cell_type {
    FLUID,
    OUTLET,
    INLET,
    NOSLIP_WALL,
    FREESLIP_WALL,
    DEFAULT
};

/// MPI coordinates struct
struct Coord {
    int x;
    int y;
};

/// MPI parameters struct
struct Params {
    int iproc = 1;
    int jproc = 1;
    int world_rank;
    int world_size;
    int neighbor_left = -1;
    int neighbor_right = -1;
    int neighbor_bottom = -1;
    int neighbor_top = -1;
    int imin;
    int imax;
    int jmin;
    int jmax;
    int size_x;
    int size_y;
};

///
struct BoundaryData{
    uint32_t neighborhood = 0;
};

/// Struct for a diagonal sparse matrix
template <typename T> struct DiagonalSparseMatrix {
    int dim;
    int num_diags = 5;
    std::vector<T> data;
    std::vector<int> offsets;
};



/**
 * @brief Class to handle logging functionality.
 *
 */
class Logger {
  public:
    /**
     * @brief Constructor for a logger object
     * 
     *Set the simulation start time at construction
     */
    Logger(){
        start_time = std::chrono::steady_clock::now();
    };

    /**
     * @brief Print text to stdout
     * 
     */
    void log(const char *str) { 
        printf("%s\n", str);
        fflush(stdout);
    }

    /**
     * @brief Print text to stderr
     * 
     */
    void log_error(const char *str) {
        fprintf(stderr, "%s\n", str);
        fflush(stderr);
    }

    /**
     * @brief Check command line argument for -log
     * 
     * Set log flag which specifies if log file is created
     * 
     * @param[in] arg command line arguments to check
     */
    void parse_flag(char *arg) {
        std::string arg_string = static_cast<std::string>(arg);
        if (arg_string == "-log") {
            _log = true;
        } else {
            std::cout << "Second argument invalid - ignored..." << std::endl;
        }
    }

    /**
     * @brief Create and write header of a log file
     * 
     * @param[in] dict name of output folder
     * @param[in] name name of the case
     * @param[in] params MPI params struct
     */
    void create_log(const std::string &dict, const std::string &name, Params &params) {
        if (_log) {
#ifdef WIN32
            char sep = '\\';
#else
            char sep = '/';
#endif
            std::string fileName = dict + sep + name + ".log";
            _log_file.open(fileName);
            _log_file << "#### " << name << " ####\n" << std::endl;
            _log_file << "Number of processes: " << params.world_size << std::endl;
            _log_file << "iproc: " << params.iproc << "   jproc: " << params.jproc << std::endl << std::endl;
        }
    }

    /**
     * @brief Print a progress bar
     * 
     * @param[in] t current simulation time
     * @param[in] t_end simulation end time
     */

    void progress_bar(double t, double t_end) {
        double progress = t / t_end;
        int pos = static_cast<int>(progress * _bar_width);
        std::cout << "Time: ";
        for (int i = 0; i < _bar_width; ++i) {
#ifndef WIN32
            if (i <= pos)
                std::cout << "\u25A0";
            else
                std::cout << "\u25A1";
#else
            if (i <= pos)
                std::cout << "=";
            else
                std::cout << "-";
#endif
        }
        std::cout << "|";
        std::cout << std::setw(5) << std::fixed << std::setprecision(2) << t;
        std::cout << " / " << t_end << "\r";
        std::cout << std::flush;
    }

    /**
     * @brief Set flag that max number of iterations were reached
     * 
     */
    void max_iter_warning() { _max_iters_reached = true; }

    /**
     * @brief Write information of current timestep to the log file
     * 
     * @param[in] timestep
     * @param[in] t simulation time
     * @param[in] dt current timestep size
     * @param[in] it pressure solver iterations
     * @param[in] max_iter
     * @param[in] res final residual
     */
    void write_log(uint32_t timestep, Real t, Real dt, int it, int max_iter, Real res) {
        if (_log) {
            _log_file << "Timestep " << timestep << ": " << std::endl;
            _log_file << "Simulation time: t = " << t << std::endl;
            _log_file << "dt: " << dt << std::endl;
            _log_file << "Final residual: " << res << "  No iterations: " << it << std::endl;
            _log_file << std::endl;
        }
    }

    /**
     * @brief log information at the end of the simulation
     * 
     * print warning about maximum iterations and runtime information
     */
    void finish() {
        auto end_time = std::chrono::steady_clock::now();
        std::chrono::duration<double> duration = end_time - start_time;

        std::cout << "\nSimulation finished in " << duration.count() << "s!" << std::endl;
        if (_max_iters_reached) {
            std::cout << "WARNING: Maximum number of iterations was reached in some timesteps" << std::endl;
        }
        if (_log) {
            _log_file << "\nTotal runtime: " << duration.count() << "s" << std::endl;
        }
    }

    // Close logFile if it was opened
    ~Logger() {
        if (_log_file.is_open()) {
            _log_file.close();
        }
    }

  private:
    // Start time for the simulation
    std::chrono::steady_clock::time_point start_time;
    // True if log file is being created
    bool _log = false;
    // Bar width of the progress bar
    const int _bar_width = 50;
    // Remember if maximum number of iterations was reached during SOR
    bool _max_iters_reached = false;
    // log File object
    std::ofstream _log_file;
};


/// Extract geometry from pgm file and create geometrical data
std::vector<std::vector<int>> parse_geometry_file(std::string filedoc, int xdim, int ydim);
/**
 * @brief Default lid driven cavity case generator
 *
 * This function creates default lid driven cavity
 * case without need for a pgm file
 */
std::vector<std::vector<int>> build_lid_driven_cavity(int xdim, int ydim);

/**
 * @brief Partition the domain for MPI processes
 * 
 * @param[in] vec complete geometry data
 * @param[in] imin local minimum i index
 * @param[in] imax local maximum i index
 * @param[in] jmin local minimum j index
 * @param[in] jmax local maximum j index
 * @param[out] local geometry data
 */
std::vector<std::vector<int>> partition(const std::vector<std::vector<int>> &vec, int imin, int imax, int jmin,
                                        int jmax);

/**
 * @brief refine the geometry by powers of 2
 * 
 * @param[in] vec geometry data
 * @param[in] refine geometry by 2 to the power of refine
 * @param[in] imax max index in x
 * @param[in] jmmax max index in y
 */
std::vector<std::vector<int>> refine_geometry(const std::vector<std::vector<int>> &vec, int refine, int &imax, int &jmax);

/// Check if cell id is Inlet id
bool is_inlet(int id);
/// Check if cell id is Outlet id
bool is_outlet(int id);
/// Check if cell id is NoSlip wall id
bool is_no_slip(int id);
/// Check if cell id is FreeSlip wall id
bool is_free_slip(int id);