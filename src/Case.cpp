#include "Case.hpp"
#include "Communication.hpp"
#include "Utilities.hpp"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
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

Case::Case(std::string file_name, int argn, char **args, Params &params) {

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
    Real KI = REAL_MAX;
    Real EPSI = REAL_MAX;
    int solver = 0;
    int refine = 0;
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
                if (var == "refine") file >> refine;
                if (var == "KI") file >> KI;
                if (var == "EPSI") file >> EPSI;
                if (var == "solver") file >> solver;
                if (var == "model") file >> _turb_model;
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

    // Check if this case uses energy equation
    if (TI != REAL_MAX) {
        _calc_temp = true;
    }

    // Prandtl number = nu / alpha
    if (pr != REAL_MAX) {
        alpha = nu / pr;
    } else if (alpha == REAL_MAX) {
        if (params.world_rank == 0 && _calc_temp) {
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

    build_domain(domain, params.size_x, params.size_y);

    _grid = Grid(_geom_name, domain, local_geometry);
    switch (_turb_model) {
    case 0: {
        std::cout << "Turbulence model: off" << std::endl;
        _field = 
            std::make_unique<Fields>(nu, dt, tau, _grid.domain().size_x, _grid.domain().size_y, 
                                     UI, VI, PI, TI, KI, EPSI, alpha, beta, GX, GY);
        break;
    }
    case 1: {
        std::cout << "Turbulence model: K-Epsilon" << std::endl;
        _field = 
            std::make_unique<TurbulenceFields>(nu, dt, tau, _grid.domain().size_x, _grid.domain().size_y, 
                                     UI, VI, PI, TI, KI, EPSI, alpha, beta, GX, GY);
        break;
    }
    default:
        break;
    }

    _discretization = Discretization(domain.dx, domain.dy, gamma);
    _max_iter = itermax;
    _tolerance = eps;

    // Construct boundaries
    if (!_grid.noslip_wall_cells().empty()) {
        _boundaries.push_back(std::make_unique<NoSlipWallBoundary>(&_grid.noslip_wall_cells(), wall_vels, wall_temps));
    }
    if (!_grid.freeslip_wall_cells().empty()) {
        _boundaries.push_back(
            std::make_unique<FreeSlipWallBoundary>(&_grid.freeslip_wall_cells(), wall_vels, wall_temps));
    }
    if (!_grid.outlet_cells().empty()) {
        _boundaries.push_back(std::make_unique<OutletBoundary>(&_grid.outlet_cells()));
    }
    if (!_grid.inlet_cells().empty()) {
        _boundaries.push_back(std::make_unique<InletBoundary>(&_grid.inlet_cells(), inlet_Us, inlet_Vs, inlet_Ts,
                                                              inlet_Ks, inlet_EPSs, DP));
    }

    switch (solver) {
    case 0: {

        _pressure_solver = std::make_unique<SOR>(omg);
        break;
    }
    case 1: {
        _pressure_solver =
            std::make_unique<PCG>(_grid.imaxb(), _grid.jmaxb(), _grid.dx(), _grid.dy(), *_field, _grid, _boundaries);
        break;
    }
    default:
        break;
    }

   
}

void Case::set_file_names(std::string file_name) {
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

/**
 * This function is the main simulation loop. In the simulation loop, following steps are required
 * - Calculate and apply boundary conditions for all the boundaries in _boundaries container
 *   using enforce_*() member function of Boundary class
 * - Calculate fluxes (F and G) using calculate_fluxes() member function of Fields class.
 *   Flux consists of diffusion and convection part, which are located in Discretization class
 * - Calculate right-hand-side of PPE using calculate_rs() member function of Fields class
 * - Iterate the pressure poisson equation until the residual becomes smaller than the desired tolerance
 *   or the maximum number of the iterations are performed using solve() member function of PressureSolver class
 * - Calculate the velocities u and v using calculate_velocities() member function of Fields class
 * - Calculate the maximal timestep size for the next iteration using calculate_dt() member function of Fields class
 * - Write vtk files using output_vtk() function
 *
 * Please note that some classes such as PressureSolver, Boundary are abstract classes which means they only provide the
 * interface. No member functions should be defined in abstract classes. You need to define functions in inherited
 * classes such as MovingWallBoundary class.
 *
 * For information about the classes and functions, you can check the header files.
 */
void Case::simulate(Params &params) {
    Real t = 0.0;
    Real dt = _field->dt();
    uint32_t timestep = 0;
    Real output_counter = 0.0;

    while (t < _t_end) {
        // Print progress bar
        if (params.world_rank == 0) logger.progress_bar(t, _t_end);

        // Select dt
        dt = _field->calculate_dt(_grid, _calc_temp);

        // Enforce velocity boundary conditions
        for (auto &boundary : _boundaries) {
            boundary->enforce_uv(*_field);
        }

        if (_calc_temp) {
            // Enforce temperature boundary conditions
            for (const auto &boundary : _boundaries) {
                boundary->enforce_t(*_field);
            }
            // Compute temperatures
            _field->calculate_temperatures(_grid);

            // Communicate temperatures
            Communication::communicate(&params, _field->t_matrix());
        }

        // Compute F & G and enforce boundary conditions
        _field->calculate_fluxes(_grid, _calc_temp);
        for (const auto &boundary : _boundaries) {
            boundary->enforce_fg(*_field);
        }
        // Communicate F and G
        Communication::communicate(&params, _field->f_matrix());
        Communication::communicate(&params, _field->g_matrix());

        // Set RHS of PPE
        _field->calculate_rs(_grid);
        // Perform pressure solve
        uint32_t it = 0;
        Real res = REAL_MAX;
        res = _pressure_solver->solve(*_field, _grid, _boundaries, params, _max_iter, _tolerance, it);

       
        // Check if max_iter was reached
        if (params.world_rank == 0 && it == _max_iter) {
            logger.max_iter_warning();
        }
        // Output current timestep information
        logger.write_log(timestep, t, dt, it, _max_iter, res);

        // Compute u^(n+1) & v^(n+1)
        _field->calculate_velocities(_grid);
        // Communicate velocities
        Communication::communicate(&params, _field->u_matrix());
        Communication::communicate(&params, _field->v_matrix());

        if(_turb_model == 1) {
            // Compute turbulent viscosity and set boundary conditions
            _field->calculate_nu_t(_grid);
            for (const auto &boundary : _boundaries) {
                boundary->enforce_nu_t(*_field);
            }
            // Communicate turbulence quantities
            Communication::communicate(&params, _field->nu_t_matrix());
            Communication::communicate(&params, _field->k_matrix());
            Communication::communicate(&params, _field->eps_matrix());
        }

        // Output u,v,p
        if (t >= output_counter * _output_freq) {
            output_vtk(timestep, params);
            output_counter++;
        }

        t += dt;
        timestep++;
        // output_vtk(timestep, params); // output every timestep for debugging
    }
    // Print Summary
    if (params.world_rank == 0) logger.finish();
    // Output u,v,p
    output_vtk(timestep, params);
}

void Case::output_vtk(int timestep, Params &params) {
    // Create a new structured grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();

    // Create grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    Real dx = _grid.dx();
    Real dy = _grid.dy();
    int i = params.world_rank % params.iproc;
    int j = params.world_rank / params.iproc;

    Real base_x = i * ((int)(global_size_x / params.iproc)) * dx + dx;
    Real base_y = j * ((int)(global_size_y / params.jproc)) * dy + dy;

    Real z = 0;
    Real y = base_y;
    for (int col = 0; col < _grid.domain().size_y + 1; col++) {
        Real x = base_x;
        for (int row = 0; row < _grid.domain().size_x + 1; row++) {
            points->InsertNextPoint(x, y, z);
            x += dx;
        }
        y += dy;
    }

    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(_grid.domain().size_x + 1, _grid.domain().size_y + 1, 1);
    structuredGrid->SetPoints(points);

// Pressure Array
#if USE_FLOATS
    typedef vtkFloatArray VTK_Array;
#else
    typedef vtkDoubleArray VTK_Array;
#endif
    VTK_Array *Pressure = VTK_Array::New();
    Pressure->SetName("pressure");
    Pressure->SetNumberOfComponents(1);

    VTK_Array *KValue = VTK_Array::New();
    KValue->SetName("kvalue");
    KValue->SetNumberOfComponents(1);

    VTK_Array *EpsValue = VTK_Array::New();
    EpsValue->SetName("epsvalue");
    EpsValue->SetNumberOfComponents(1);

    VTK_Array *TurbViscosity = VTK_Array::New();
    TurbViscosity->SetName("nu_t");
    TurbViscosity->SetNumberOfComponents(1);
    // Velocity Array
    VTK_Array *Velocity = VTK_Array::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);

    // Print pressure and place ghost cells
    for (int j = 1; j < _grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _grid.domain().size_x + 1; i++) {
            Real pressure = _field->p(i, j);
            Pressure->InsertNextTuple(&pressure);
            Real kval = _field->k(i, j);
            KValue->InsertNextTuple(&kval);
            Real epsval = _field->eps(i, j);
            EpsValue->InsertNextTuple(&epsval);
            Real nu_t = _field->nu_t(i, j);
            TurbViscosity->InsertNextTuple(&nu_t);
            // Insert blank cells at obstacles
            if (_grid.cell(i, j).type() != cell_type::FLUID) {
                structuredGrid->BlankCell((j - 1) * _grid.domain().domain_size_x + (i - 1));
            }
        }
    }

    // Temp Velocity
    float vel[3];
    vel[2] = 0; // Set z component to 0

    // Print Velocity from bottom to top
    for (int j = 0; j < _grid.domain().size_y + 1; j++) {
        for (int i = 0; i < _grid.domain().size_x + 1; i++) {
            vel[0] = (_field->u(i, j) + _field->u(i, j + 1)) * 0.5;
            vel[1] = (_field->v(i, j) + _field->v(i + 1, j)) * 0.5;
            Velocity->InsertNextTuple(vel);
        }
    }

    // Add Pressure to Structured Grid
    structuredGrid->GetCellData()->AddArray(Pressure);

    structuredGrid->GetCellData()->AddArray(KValue);
    structuredGrid->GetCellData()->AddArray(EpsValue);
    structuredGrid->GetCellData()->AddArray(TurbViscosity);

    // Add Velocity to Structured Grid
    structuredGrid->GetPointData()->AddArray(Velocity);

    // Add Temperature to Structured Grid
    if (_calc_temp) {
        VTK_Array *Temperature = VTK_Array::New();
        Temperature->SetName("temperature");
        Temperature->SetNumberOfComponents(1);

        for (int j = 1; j < _grid.domain().size_y + 1; j++) {
            for (int i = 1; i < _grid.domain().size_x + 1; i++) {

                Real temperature = _field->t(i, j);
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

void Case::build_domain(Domain &domain, int imax_domain, int jmax_domain) {
    domain.imin = 0;
    domain.jmin = 0;
    domain.imax = imax_domain + 2;
    domain.jmax = jmax_domain + 2;
    domain.size_x = imax_domain;
    domain.size_y = jmax_domain;
}
