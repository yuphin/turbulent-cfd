#include "Simulation.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>
namespace filesystem = std::filesystem;
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkTuple.h>

Simulation::Simulation(std::string file_name, int argn, char **args, Params &params) {

    // Set up logging functionality
    if (argn > 2) {
        logger.parse_flag(args[2]);
    }

    // Read input parameters
    const int MAX_LINE_LENGTH = 1024;
    std::ifstream file(file_name);
    Real nu = REAL_MAX;    /* viscosity   */
    Real UI;               /* velocity x-direction */
    Real VI;               /* velocity y-direction */
    Real PI;               /* pressure */
    Real GX;               /* gravitation x-direction */
    Real GY;               /* gravitation y-direction */
    Real xlength;          /* length of the domain x-dir.*/
    Real ylength;          /* length of the domain y-dir.*/
    Real dt;               /* time step */
    int imax;              /* number of cells x-direction*/
    int jmax;              /* number of cells y-direction*/
    Real gamma;            /* uppwind differencing factor*/
    Real omg;              /* relaxation factor */
    Real tau;              /* safety factor for time step*/
    int itermax;           /* max. number of iterations for pressure per time step */
    Real eps;              /* accuracy bound for pressure*/
    Real re = REAL_MAX;    /* Reynolds number */
    Real pr = REAL_MAX;    /* Prandtl number */
    Real beta;             /* thermal expansion coefficient */
    Real TI = REAL_MAX;    /* Temperature */
    Real alpha = REAL_MAX; /* Thermal diffusivity */
    Real DP = REAL_MAX;    /* Pressure differential between the two ends */
    int solver_type_int = 0;
    int simulation_type_int = 0;
    int preconditioner = -1;
    SolverType solver_type;
    Real KI = REAL_MAX;
    Real EPSI = REAL_MAX;
    int refine = 0;        /* Refinement based on power of 2 */
    int turb_model = 0;
    std::unordered_map<int, Real> wall_temps;
    std::unordered_map<int, Real> wall_vels;
    std::unordered_map<int, Real> inlet_Us;
    std::unordered_map<int, Real> inlet_Vs;
    std::unordered_map<int, Real> inlet_Ts;
    std::unordered_map<int, Real> inlet_Ks;
    std::unordered_map<int, Real> inlet_EPSs;
    if (file.is_open()) {
        std::string var;
        while (!file.eof() && file.good()) {
            file >> var;
            if (var[0] == '#') { /* ignore comment line*/
                file.ignore(MAX_LINE_LENGTH, '\n');
            } else {
                if (var == "geo_file") file >> _geom_name;
                if (var == "xlength") file >> xlength;
                if (var == "ylength") file >> ylength;
                if (var == "nu") file >> nu;
                if (var == "t_end") file >> _t_end;
                if (var == "dt") file >> dt;
                if (var == "omg") file >> omg;
                if (var == "eps") file >> eps;
                if (var == "tau") file >> tau;
                if (var == "gamma") file >> gamma;
                if (var == "dt_value") file >> _output_freq;
                if (var == "UI") file >> UI;
                if (var == "VI") file >> VI;
                if (var == "GX") file >> GX;
                if (var == "GY") file >> GY;
                if (var == "PI") file >> PI;
                if (var == "itermax") file >> itermax;
                if (var == "imax") file >> imax;
                if (var == "jmax") file >> jmax;
                if (var == "Re") file >> re;
                if (var == "Pr") file >> pr;
                if (var == "beta") file >> beta;
                if (var == "TI") file >> TI;
                if (var == "alpha") file >> alpha;
                if (var == "DELTA_P") file >> DP;
                if (var == "iproc") file >> params.iproc;
                if (var == "jproc") file >> params.jproc;
                if (var == "solver") file >> solver_type_int;
                if (var == "simulation") file >> simulation_type_int;
                if (var == "preconditioner") file >> preconditioner;
                if (var == "refine") file >> refine;
                if (var == "model") file >> turb_model;
                if (var == "KI") file >> KI;
                if (var == "EPSI") file >> EPSI;
                if (!var.compare(0, 10, "wall_temp_")) {
                    Real temp;
                    file >> temp;
                    if (temp == -1.0) {
                    }
                    wall_temps.insert({std::stoi(var.substr(10)), temp});
                }
                if (!var.compare(0, 9, "wall_vel_")) {
                    Real vel;
                    file >> vel;
                    wall_vels.insert({std::stoi(var.substr(9)), vel});
                }
                if (!var.compare(0, 4, "UIN_")) {
                    Real u;
                    file >> u;
                    inlet_Us.insert({std::stoi(var.substr(4)), u});
                }
                if (!var.compare(0, 4, "VIN_")) {
                    Real v;
                    file >> v;
                    inlet_Vs.insert({std::stoi(var.substr(4)), v});
                }
                if (!var.compare(0, 4, "TIN_")) {
                    Real t;
                    file >> t;
                    inlet_Ts.insert({std::stoi(var.substr(4)), t});
                }
                if (!var.compare(0, 4, "KIN_")) {
                    Real k;
                    file >> k;
                    inlet_Ks.insert({std::stoi(var.substr(4)), k});
                }
                if (!var.compare(0, 6, "EPSIN_")) {
                    Real eps;
                    file >> eps;
                    inlet_EPSs.insert({std::stoi(var.substr(6)), eps});
                }
            }
        }
    }
    file.close();
    if (params.iproc * params.jproc != params.world_size) {
        if (params.world_rank == 0)
            std::cout << "ERROR: Number of MPI processes doesn't match iproc * jproc! \nAborting... " << std::endl;
        Communication::finalize();
        std::exit(0);
    }
    // We assume Reynolds number = 1 / nu for now
    if (re != REAL_MAX && nu == REAL_MAX) {
        nu = 1 / re;
    } else if (re == REAL_MAX && nu == REAL_MAX) {
        if (params.world_rank == 0) {
            logger.log_error("Viscosity and Reynolds number not specified, defaulting viscosity to 0");
        }
        nu = 0.0;
    }
    if (solver_type_int == (int)SolverType::SOR) {
        solver_type = SolverType::SOR;
    } else if (solver_type_int == (int)SolverType::PCG) {
        solver_type = SolverType::PCG;
    }

    if (simulation_type_int == 0) {
        _solver = std::make_unique<CPUSolver>();
        if (params.world_rank == 0) {
            std::cout << "Simulation: CPU\n";
        }
    }
#ifdef USE_CUDA
    else if (simulation_type_int == 1) {
        _solver = std::make_unique<CudaSolver>();
        std::cout << "Simulation: CUDA\n";
    }
#endif
#ifdef USE_VULKAN
    else if (simulation_type_int == 2) {
        _solver = std::make_unique<VulkanSolver>();
        std::cout << "Simulation: Vulkan\n";
    }
#endif

    if (params.world_rank == 0) {
        if (sizeof(Real) == 4) {
            std::cout << "Precision: Single\n";
        } else {
            std::cout << "Precision: Double\n";
        }        
    }

    // Prandtl number = nu / alpha
    if (pr != REAL_MAX) {
        alpha = nu / pr;
    } else if (alpha == REAL_MAX) {
        if (params.world_rank == 0) {
            logger.log_error("Prandtl number, alpha or beta are not set, defaulting to 0");
        }
        alpha = 0.0;
        beta = 0.0;
    }

    if (_geom_name.compare("NONE") == 0) {
        wall_vels.insert(std::pair<int, Real>(LidDrivenCavity::moving_wall_id, LidDrivenCavity::wall_velocity));
    }

    // Set file names for geometry file and output directory
    set_file_names(file_name);
    // Create log file in output dir
    logger.create_log(_dict_name, _case_name, params);

    std::vector<std::vector<int>> global_geometry;
    if (_geom_name.compare("NONE")) {

        global_geometry = parse_geometry_file(_geom_name, imax, jmax);
    } else {
        global_geometry = build_lid_driven_cavity(imax, jmax);
    }

    global_geometry = refine_geometry(global_geometry, refine, imax, jmax);
    global_size_x = imax;
    global_size_y = jmax;
    Communication::init_params(&params, imax, jmax);

    auto local_geometry = partition(global_geometry, params.imin, params.imax, params.jmin, params.jmax);
    // Build up the domain
    Domain domain;
    domain.dx = xlength / (Real)imax;
    domain.dy = ylength / (Real)jmax;
    domain.domain_size_x = params.size_x;
    domain.domain_size_y = params.size_y;
    domain.x_length = xlength;
    domain.y_length = ylength;

    domain.imin = 0;
    domain.jmin = 0;
    domain.imax = params.size_x + 2;
    domain.jmax = params.size_y + 2;
    domain.size_x = params.size_x;
    domain.size_y = params.size_y;
    domain.total_size = domain.imax * domain.jmax;

    _solver->_grid = Grid(domain, local_geometry);
    _solver->_discretization = Discretization(domain.dx, domain.dy, gamma);
    _solver->_max_iter = itermax;
    _solver->_tolerance = eps;
    _solver->params = params;
    _solver->_field = Fields(nu, dt, tau, domain.size_x, domain.size_y, UI, VI, PI, TI, KI, EPSI, alpha, beta, GX, GY);
    _solver->solver_type = solver_type;
    _solver->_preconditioner = preconditioner;
    _solver->_omega = omg;
    _solver->_turb_model = turb_model;
    // Check if this case uses energy equation
    if (TI != REAL_MAX) {
        _solver->_field.calc_temp = true;
    }
    if (params.world_rank == 0) {
        if (turb_model == 0) {
            std::cout << "Turbulence model: off\n";
        } else if (turb_model == 1) {
            std::cout << "Turbulence model: K-Epsilon\n";
        } else if (turb_model == 2) {
            std::cout << "Turbulence model: K-Omega\n";
        } else if (turb_model == 3) {
            std::cout << "Turbulence model: K-Omega SST\n";
        }        
    }

    // Construct boundaries
    if (!_solver->_grid.outlet_cells().empty()) {
        _solver->_boundaries.push_back(std::make_unique<OutletBoundary>(&_solver->_grid.outlet_cells()));
    }
    if (!_solver->_grid.inlet_cells().empty()) {
        _solver->_boundaries.push_back(std::make_unique<InletBoundary>(&_solver->_grid.inlet_cells(), inlet_Us,
                                                                       inlet_Vs, inlet_Ts, inlet_Ks, inlet_EPSs, DP));
    }
    if (!_solver->_grid.freeslip_wall_cells().empty()) {
        _solver->_boundaries.push_back(
            std::make_unique<FreeSlipWallBoundary>(&_solver->_grid.freeslip_wall_cells(), wall_vels, wall_temps));
    }

    if (!_solver->_grid.noslip_wall_cells().empty()) {
        _solver->_boundaries.push_back(
            std::make_unique<NoSlipWallBoundary>(&_solver->_grid.noslip_wall_cells(), wall_vels, wall_temps));
    }

    // TODO
    for (auto &p : inlet_EPSs) {
        _solver->_EPSIN = p.second;
    }
    for (auto &p : inlet_Ks) {
        _solver->_KIN = p.second;
    }
}

void Simulation::set_file_names(std::string file_name) {
    std::string temp_dir;
    bool case_name_flag = true;
    bool prefix_flag = false;

    for (int i = static_cast<int>(file_name.size()) - 1; i > -1; --i) {
        if (file_name[i] == '/') {
            case_name_flag = false;
            prefix_flag = true;
        }
        if (case_name_flag) {
            _case_name.push_back(file_name[i]);
        }
        if (prefix_flag) {
            _prefix.push_back(file_name[i]);
        }
    }

    for (int i = static_cast<int>(file_name.size() - _case_name.size() - 1); i > -1; --i) {
        temp_dir.push_back(file_name[i]);
    }

    std::reverse(_case_name.begin(), _case_name.end());
    std::reverse(_prefix.begin(), _prefix.end());
    std::reverse(temp_dir.begin(), temp_dir.end());

    _case_name.erase(_case_name.size() - 4);
    _dict_name = temp_dir;
    _dict_name.append(_case_name);
    _dict_name.append("_Output");

    if (_geom_name.compare("NONE") != 0) {
        _geom_name = _prefix + _geom_name;
    }

    // Create output directory
    filesystem::path folder(_dict_name);
    try {
        filesystem::remove_all(folder);
        filesystem::create_directory(folder);
    } catch (const std::exception &) {
        std::cerr << "Output directory could not be created." << std::endl;
        std::cerr << "Make sure that you have write permissions to the "
                     "corresponding location"
                  << std::endl;
    }
}

void Simulation::simulate(Params &params) {
    Real t = 0.0;
    Real dt;
    uint32_t timestep = 0;
    Real output_counter = 0.0;
    _solver->initialize();
    while (t < _t_end) {
        // Print progress bar
        if (params.world_rank == 0) logger.progress_bar(t, _t_end);

        _solver->solve_pre_pressure(dt);
        uint32_t it;
        Real res;
        _solver->solve_pressure(res, it);
        
        if (params.world_rank == 0) {
            std::cout << "Iter count: " << it << " ";
        }
        
        // Check if max_iter was reached
        if (params.world_rank == 0 && it == _solver->_max_iter) {
            logger.max_iter_warning();
        }
        // Output current timestep information
        logger.write_log(timestep, t, dt, it, _solver->_max_iter, res);
        _solver->solve_post_pressure();
        // Output u,v,p
        if (t >= output_counter * _output_freq) {
            output_vtk(timestep, params);
            output_counter++;
        }

        t += dt;
        timestep++;
        // output_vtk(timestep, params); // output every timestep for debugging
    }
    _solver->solve_post_pressure();
}

void Simulation::output_vtk(int timestep, Params &params) {
    // Create a new structured grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();

    // Create grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    Real dx = _solver->_grid.dx();
    Real dy = _solver->_grid.dy();
    int i = params.world_rank % params.iproc;
    int j = params.world_rank / params.iproc;

    Real base_x = i * ((int)(global_size_x / params.iproc)) * dx + dx;
    Real base_y = j * ((int)(global_size_y / params.jproc)) * dy + dy;

    Real z = 0;
    Real y = base_y;
    for (int col = 0; col < _solver->_grid.domain().size_y + 1; col++) {
        Real x = base_x;
        for (int row = 0; row < _solver->_grid.domain().size_x + 1; row++) {
            points->InsertNextPoint(x, y, z);
            x += dx;
        }
        y += dy;
    }

    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(_solver->_grid.domain().size_x + 1, _solver->_grid.domain().size_y + 1, 1);
    structuredGrid->SetPoints(points);

// Pressure Array
#if USE_FLOATS
    typedef vtkFloatArray VTK_Array;
#else
    typedef vtkDoubleArray VTK_Array;
#endif
    vtkSmartPointer<VTK_Array> Pressure = vtkSmartPointer<VTK_Array>::New();
    Pressure->SetName("pressure");
    Pressure->SetNumberOfComponents(1);
    vtkSmartPointer<VTK_Array> KValue;
    vtkSmartPointer<VTK_Array> EpsValue;
    vtkSmartPointer<VTK_Array> TurbViscosity;
    if (_solver->_turb_model != 0) {
        KValue = vtkSmartPointer<VTK_Array>::New();
        KValue->SetName("kvalue");
        KValue->SetNumberOfComponents(1);

        EpsValue = vtkSmartPointer<VTK_Array>::New();
        EpsValue->SetName("epsvalue");
        EpsValue->SetNumberOfComponents(1);

        TurbViscosity = vtkSmartPointer<VTK_Array>::New();
        TurbViscosity->SetName("nu_t");
        TurbViscosity->SetNumberOfComponents(1);
    }

    // Velocity Array
    vtkSmartPointer<VTK_Array> Velocity = vtkSmartPointer<VTK_Array>::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);

    // Print pressure and place ghost cells
    for (int j = 1; j < _solver->_grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _solver->_grid.domain().size_x + 1; i++) {
            Real pressure = _solver->_field.p(i, j);
            Pressure->InsertNextTuple(&pressure);
            if (_solver->_turb_model != 0) {
                Real kval = _solver->_field.k(i, j);
                KValue->InsertNextTuple(&kval);
                Real epsval = _solver->_field.eps(i, j);
                EpsValue->InsertNextTuple(&epsval);
                Real nu_t = _solver->_field.nu_t(i, j);
                TurbViscosity->InsertNextTuple(&nu_t);
            }
            // Insert blank cells at obstacles
            if (_solver->_grid.cell(i, j).type() != cell_type::FLUID) {
                structuredGrid->BlankCell((j - 1) * _solver->_grid.domain().domain_size_x + (i - 1));
            }
        }
    }

    // Temp Velocity
    float vel[3];
    vel[2] = 0; // Set z component to 0

    // Print Velocity from bottom to top
    for (int j = 0; j < _solver->_grid.domain().size_y + 1; j++) {
        for (int i = 0; i < _solver->_grid.domain().size_x + 1; i++) {
            vel[0] = (_solver->_field.u(i, j) + _solver->_field.u(i, j + 1)) * 0.5;
            vel[1] = (_solver->_field.v(i, j) + _solver->_field.v(i + 1, j)) * 0.5;
            Velocity->InsertNextTuple(vel);
        }
    }

    // Add Pressure to Structured Grid
    structuredGrid->GetCellData()->AddArray(Pressure);
    if (_solver->_turb_model != 0) {
        structuredGrid->GetCellData()->AddArray(KValue);
        structuredGrid->GetCellData()->AddArray(EpsValue);
        structuredGrid->GetCellData()->AddArray(TurbViscosity);
    }
    // Add Velocity to Structured Grid
    structuredGrid->GetPointData()->AddArray(Velocity);

    // Add Temperature to Structured Grid
    if (_solver->_field.calc_temp) {
        vtkSmartPointer<VTK_Array> Temperature = VTK_Array::New();
        Temperature->SetName("temperature");
        Temperature->SetNumberOfComponents(1);

        for (int j = 1; j < _solver->_grid.domain().size_y + 1; j++) {
            for (int i = 1; i < _solver->_grid.domain().size_x + 1; i++) {

                Real temperature = _solver->_field.t(i, j);
                Temperature->InsertNextTuple(&temperature);
            }
        }
        structuredGrid->GetCellData()->AddArray(Temperature);
    }

    // Write Grid
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();

    // Create Filename
    std::string outputname = _dict_name + '/' + _case_name + "_" + std::to_string(params.world_rank) + "." +
                             std::to_string(timestep) + ".vtk";

    writer->SetFileName(outputname.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
}
