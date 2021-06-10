#pragma once

#include <cfloat>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#define USE_FLOATS 1

#if USE_FLOATS
const float REAL_MAX = FLT_MAX;
typedef float Real;
#else
const double REAL_MAX = DBL_MAX;
typedef double Real;
#endif

struct BoundaryData {
    uint32_t neighborhood = 0;
};


template <typename T>
struct DiagonalSparseMatrix {
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
    Logger(){};

    // Check flag for creating a log file
    void parse_flag(char *arg) {
        std::string arg_string = static_cast<std::string>(arg);
        if (arg_string == "-log") {
            _log = true;
        } else {
            std::cout << "Second argument invalid - ignored..." << std::endl;
        }
    }

    // Prepare the log file
    void create_log(const std::string &dict, const std::string &name) {
        if (_log) {
#ifdef WIN32
            char sep = '\\';
#else
            char sep = '/';
#endif
            std::string fileName = dict + sep + name + ".log";
            _log_file.open(fileName);
            _log_file << "#### " << name << " ####\n" << std::endl;
        }
    }

    // Print a progress bar
    void progress_bar(double t, double t_end) {
        double progress = t / t_end;
        int pos = static_cast<int>(progress * _bar_width);
        std::cout << " Time: |";
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

    // Set _maxItersReached
    void max_iter_warning() { _max_iters_reached = true; }

    // Write information of current timestep into logFile
    void write_log(uint32_t timestep, Real t, int it, int max_iter, Real res) {
        if (_log) {
            _log_file << "Timestep " << timestep << ": " << std::endl;
            _log_file << "Simulation time: t = " << t << std::endl;
            _log_file << "SOR final residual: " << res << "  No iterations: " << it << std::endl;
            _log_file << std::endl;
        }
    }

    // Final output at the end of the simulation
    void finish() {
        std::cout << "\nSimulation finished!" << std::endl;
        if (_max_iters_reached) {
            std::cout << "WARNING: Maximum number of iterations was reached in some timesteps" << std::endl;
        }
    }

    // Close logFile if it was opened
    ~Logger() {
        if (_log_file.is_open()) {
            _log_file.close();
        }
    }

  private:
    // True if log file is being created
    bool _log = false;
    // Bar width of the progress bar
    const int _bar_width = 50;
    // Remember if maximum number of iterations was reached during SOR
    bool _max_iters_reached = false;
    // log File object
    std::ofstream _log_file;
};