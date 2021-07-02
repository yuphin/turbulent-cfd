#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Boundary.hpp"
#include "Discretization.hpp"
#include "Domain.hpp"
#include "Fields.hpp"
#include "Grid.hpp"
#include "PressureSolver.hpp"
#include "Solver.hpp"
#include "Utilities.hpp"
#include "Communication.hpp"
#include "Utilities.hpp"
#include "CPUSolver.hpp"
#include "VulkanSolver.hpp"

/**
 * @brief Class to hold and orchestrate the simulation flow.
 *
 */
class Simulation {
  public:
    /**
     * @brief Parallel constructor for the Simulation.
     *
     * Reads input file, creates Fields, Grid, Boundary, Solver and sets
     * Discretization parameters Creates output directory
     *
     * @param[in] Input file name
     */
    Simulation(std::string file_name, int argn, char **args, Params &params);

    /**
     * @brief Main function to simulate the flow until the end time.
     *
     * Calculates the fluxes
     * Calculates the right hand side
     * Solves pressure
     * Calculates velocities
     * Outputs the solution files
     */
    void simulate(Params &params);

  private:
    /// Plain case name without paths
    std::string _case_name;
    /// Output directiory name
    std::string _dict_name;
    /// Geometry file name
    std::string _geom_name{"NONE"};
    /// Relative input file path
    std::string _prefix;

    /// Simulation time
    Real _t_end;
    /// Solution file outputting frequency
    Real _output_freq;

    std::unique_ptr<Solver> _solver;

    friend class Logger;
    Logger logger = Logger();

    // global grid sizes without boundaries, needed for vtk output
    Real global_size_x;
    Real global_size_y;

    /**
     * @brief Creating file names from given input data file
     *
     * Extracts path of the case file and creates code-readable file names
     * for outputting directory and geometry file.
     *
     * @param[in] input data file
     */
    void set_file_names(std::string file_name);

    /**
     * @brief Solution file outputter
     *
     * Outputs the solution files in .vtk format. Ghost cells are excluded.
     * Pressure is cell variable while velocity is point variable while being
     * interpolated to the cell faces
     *
     * @param[in] Timestep of the solution
     */
    void output_vtk(int t, Params &params);

    void build_domain(Domain &domain, int imax_domain, int jmax_domain);
};
