#include "Case.hpp"
#include "Enums.hpp"
#include "GPUSimulation.hpp"
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

Case::Case(std::string file_name, int argn, char **args) {

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
    std::unordered_map<int, Real> wall_temps;
    std::unordered_map<int, Real> wall_vels;
    std::unordered_map<int, Real> inlet_Us;
    std::unordered_map<int, Real> inlet_Vs;
    std::unordered_map<int, Real> inlet_Ts;

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
            }
        }
    }
    file.close();

    // We assume Reynolds number = 1 / nu for now
    if (re != REAL_MAX && nu == REAL_MAX) {
        nu = 1 / re;
    } else if (re == REAL_MAX && nu == REAL_MAX) {
        std::cerr << "Viscosity and Reynolds number not specified, defaulting viscosity to 0\n";
        nu = 0.0;
    }

    // Prandtl number = nu / alpha
    if (pr != REAL_MAX) {
        alpha = nu / pr;
    } else if (alpha == REAL_MAX) {
        std::cerr << "Prandtl number, alpha or beta are not set, defaulting to 0\n";
        alpha = 0.0;
        beta = 0.0;
    }

    if (TI != REAL_MAX) {
        _calc_temp = true;
    }

    if (_geom_name.compare("NONE") == 0) {
        wall_vels.insert(std::pair<int, Real>(LidDrivenCavity::moving_wall_id, LidDrivenCavity::wall_velocity));
    }

    // Set file names for geometry file and output directory
    set_file_names(file_name);
    // Create log file in output dir
    logger.create_log(_dict_name, _case_name);

    // Build up the domain
    Domain domain;
    domain.dx = xlength / (Real)imax;
    domain.dy = ylength / (Real)jmax;
    domain.domain_size_x = imax;
    domain.domain_size_y = jmax;

    build_domain(domain, imax, jmax);

    _grid = Grid(_geom_name, domain);
    _field = Fields(nu, dt, tau, _grid.domain().size_x, _grid.domain().size_y, UI, VI, PI, TI, alpha, beta, GX, GY);

    _discretization = Discretization(domain.dx, domain.dy, gamma);
    _pressure_solver = std::make_unique<SOR>(omg);
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
        _boundaries.push_back(std::make_unique<InletBoundary>(&_grid.inlet_cells(), inlet_Us, inlet_Vs, inlet_Ts, DP));
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
void Case::simulate() {
    Real t = 0.0;
    Real dt = _field.dt();
    uint32_t timestep = 0;
    Real output_counter = 0.0;
    GPUSimulation simulation;
    simulation.init();
    simulation.create_descriptor_set_layout(
        {descriptor_set_layout_binding(VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, VK_SHADER_STAGE_COMPUTE_BIT, 0),
         descriptor_set_layout_binding(VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, VK_SHADER_STAGE_COMPUTE_BIT, 1),
         descriptor_set_layout_binding(VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, VK_SHADER_STAGE_COMPUTE_BIT, 2),
         descriptor_set_layout_binding(VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, VK_SHADER_STAGE_COMPUTE_BIT, 3),
         descriptor_set_layout_binding(VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, VK_SHADER_STAGE_COMPUTE_BIT, 4),
         descriptor_set_layout_binding(VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, VK_SHADER_STAGE_COMPUTE_BIT, 5),
         descriptor_set_layout_binding(VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, VK_SHADER_STAGE_COMPUTE_BIT, 6),
         descriptor_set_layout_binding(VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, VK_SHADER_STAGE_COMPUTE_BIT, 7),
         descriptor_set_layout_binding(VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, VK_SHADER_STAGE_COMPUTE_BIT, 8),
         descriptor_set_layout_binding(VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, VK_SHADER_STAGE_COMPUTE_BIT, 9)});
    simulation.create_command_pool();
    simulation.create_fences();
    Pipeline discretization_pipeline = simulation.create_compute_pipeline("src/shaders/discretization.comp.spv");
    Pipeline rs_pipeline = simulation.create_compute_pipeline("src/shaders/calc_rs.comp.spv");
    Pipeline vel_pipeline = simulation.create_compute_pipeline("src/shaders/calc_vel.comp.spv");
    Pipeline p_pipeline_red = simulation.create_compute_pipeline("src/shaders/calc_p_redblack_gs.comp.spv", 0);
    Pipeline p_pipeline_black = simulation.create_compute_pipeline("src/shaders/calc_p_redblack_gs.comp.spv", 1);
    Pipeline residual_pipeline = simulation.create_compute_pipeline("src/shaders/calc_r.comp.spv");
    Pipeline p_boundary_pipeline = simulation.create_compute_pipeline("src/shaders/boundary_p.comp.spv");

    std::vector<Real> is_fluid(_grid.imaxb() * _grid.jmaxb(), 0);
    for (const auto &current_cell : _grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();
        is_fluid[_grid.imaxb() * j + i] = 1;
    }
    // Preprocess
    std::vector<BoundaryData> boundaries(_grid.imaxb() * _grid.jmaxb());
    for (const auto &boundary : _boundaries) {
        BoundaryData data;
        uint32_t type = boundary->get_type();
        auto cells = boundary->_cells;
        data.neighborhood |= type << 8;
        for (auto &cell : *cells) {
            int i = cell->i();
            int j = cell->j();
            data.neighborhood |= (-1 & 0xFF);
            // data.idx = _grid.imaxb() * j + i;
            if (cell->is_border(border_position::RIGHT)) {
                data.neighborhood = (data.neighborhood & (-1 << 8)) | 0;
            }
            if (cell->is_border(border_position::LEFT)) {
                data.neighborhood = (data.neighborhood & (-1 << 8)) | 1;
            }
            if (cell->is_border(border_position::TOP)) {
                data.neighborhood = (data.neighborhood & (-1 << 8)) | 2;
            }
            if (cell->is_border(border_position::BOTTOM)) {
                data.neighborhood = (data.neighborhood & (-1 << 8)) | 3;
            }
            if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::TOP)) {
                data.neighborhood = (data.neighborhood & (-1 << 8)) | 4;
            }
            if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::BOTTOM)) {
                data.neighborhood = (data.neighborhood & (-1 << 8)) | 5;
            }
            if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::TOP)) {
                data.neighborhood = (data.neighborhood & (-1 << 8)) | 6;
            }
            if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::BOTTOM)) {
                data.neighborhood = (data.neighborhood & (-1 << 8)) | 7;
            }
            boundaries[j * _grid.imaxb() + i] = data;
        }
    }
    UBOData ubo_data;
    Buffer t_buffer;
    Buffer cell_buffer;
    Buffer neighborhood_buffer;
    Buffer u_buffer;
    Buffer v_buffer;
    Buffer f_buffer;
    Buffer g_buffer;
    Buffer rs_buffer;
    Buffer p_buffer;
    Buffer residual_buffer;
    Buffer scratch_buffer;
    t_buffer.create(&simulation.context, VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT,
                    VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                    VK_SHARING_MODE_EXCLUSIVE, sizeof(UBOData));
    cell_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                       VK_SHARING_MODE_EXCLUSIVE, sizeof(Real) * _field.f_matrix().size(), is_fluid.data(), true);

    neighborhood_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                               VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                               VK_SHARING_MODE_EXCLUSIVE, boundaries.size() * sizeof(BoundaryData), boundaries.data());
    u_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                    VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                    VK_SHARING_MODE_EXCLUSIVE, _field.u_matrix().size() * sizeof(Real), _field._U._container.data());
    v_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                    VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                    VK_SHARING_MODE_EXCLUSIVE, _field.v_matrix().size() * sizeof(Real), _field._V._container.data());
    f_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                    VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                    VK_SHARING_MODE_EXCLUSIVE, _field.f_matrix().size() * sizeof(Real), _field._F._container.data());
    g_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                    VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                    VK_SHARING_MODE_EXCLUSIVE, _field.g_matrix().size() * sizeof(Real), _field._G._container.data());

    residual_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                           VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                           VK_SHARING_MODE_EXCLUSIVE, 50 * sizeof(Real));
    rs_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
                     VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE, _field._RS.size() * sizeof(Real),
                     _field._RS._container.data(), true);
    p_buffer.create(&simulation.context,
                    VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT |
                        VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                    VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE, _field._P.size() * sizeof(Real));
    scratch_buffer.create(&simulation.context, VK_BUFFER_USAGE_TRANSFER_DST_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
                          VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                          VK_SHARING_MODE_EXCLUSIVE, _field._P.size() * sizeof(Real));
    simulation.push_descriptors({{neighborhood_buffer, 9}});
    while (t < _t_end) {

        // Print progress bar
        logger.progress_bar(t, _t_end);
        // Select dt
        dt = _field.calculate_dt(_grid, _calc_temp);
        // Enforce velocity boundary conditions
        for (auto &boundary : _boundaries) {
            boundary->enforce_uv(_field);
        }
        u_buffer.upload(_field.u_matrix().size() * sizeof(Real), _field._U._container.data());
        v_buffer.upload(_field.v_matrix().size() * sizeof(Real), _field._V._container.data());
        ubo_data.dt = dt;
        memcpy(t_buffer.data, &ubo_data, sizeof(UBOData));

        simulation.push_descriptors({{u_buffer, 0},
                                     {v_buffer, 1},
                                     {f_buffer, 2},
                                     {g_buffer, 3},
                                     {t_buffer, 4, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER},
                                     {cell_buffer, 5},
                                     {rs_buffer, 6}});

        if (_calc_temp) {
            // Enforce temperature boundary conditions
            for (const auto &boundary : _boundaries) {
                boundary->enforce_t(_field);
            }
            // Compute temperatures
            _field.calculate_temperatures(_grid);
        }

        // Compute F & G and enforce boundary conditions
        //_field.calculate_fluxes(_grid, _calc_temp);
        simulation.begin_end_record_command_buffer(discretization_pipeline);
        simulation.run_command_buffer();
        _field.f_matrix()._container.assign((Real *)f_buffer.data, (Real *)f_buffer.data + _field.f_matrix().size());
        _field.g_matrix()._container.assign((Real *)g_buffer.data, (Real *)g_buffer.data + _field.g_matrix().size());
        for (const auto &boundary : _boundaries) {
            boundary->enforce_fg(_field);
        }
        // Set RHS of PPE
        //_field.calculate_rs(_grid);
        f_buffer.upload(_field.f_matrix().size() * sizeof(Real), _field._F._container.data());
        g_buffer.upload(_field.g_matrix().size() * sizeof(Real), _field._G._container.data());

        simulation.push_descriptors({{f_buffer, 2}, {g_buffer, 3}});
        simulation.begin_end_record_command_buffer(rs_pipeline);
        simulation.run_command_buffer();
        // Perform pressure solve
        uint32_t it = 0;
        Real res = REAL_MAX;
        //_max_iter = 100;
        rs_buffer.copy(scratch_buffer);
        _field._RS._container.assign((Real *)scratch_buffer.data, (Real *)scratch_buffer.data + _field._RS.size());
        // Do some sequential Gauss Seidel Iterations
        while (it < 100 && res > _tolerance) {
            res = _pressure_solver->solve(_field, _grid, _boundaries);
            // Enforce boundary conditions
            for (const auto &boundary : _boundaries) {
                boundary->enforce_p(_field);
            }
            it++;
        }
        it = 0;
        p_buffer.upload(_field.p_matrix().size() * sizeof(Real), _field._P._container.data(), &scratch_buffer);
        simulation.push_descriptors({{p_buffer, 7}, {residual_buffer, 8}});
        // Follow it by Red-Black Gauss Seidel
        while (it < _max_iter && res > _tolerance) {
            simulation.begin_recording();
            simulation.record_command_buffer(p_pipeline_red);
            simulation.record_command_buffer(p_pipeline_black);
            simulation.record_command_buffer(residual_pipeline);
            simulation.record_command_buffer(p_boundary_pipeline);
            simulation.end_recording();
            simulation.run_command_buffer();
            it++;
        }
        // You can uncomment this for debugging
        /* p_buffer.copy(scratch_buffer);
         _field._P._container.assign((Real *)scratch_buffer.data, (Real *)scratch_buffer.data + _field._P.size());*/

        // Check if max_iter was reached
        if (it == _max_iter) {
            // std::cout << it << " " << res;
            logger.max_iter_warning();
        }
        // Output current timestep information
        logger.write_log(timestep, t, it, _max_iter, res);

        // Compute u^(n+1) & v^(n+1)
        //_field.calculate_velocities(_grid);

        simulation.begin_end_record_command_buffer(vel_pipeline);
        simulation.run_command_buffer();
        _field._U._container.assign((Real *)u_buffer.data, (Real *)u_buffer.data + _field._U.size());
        _field._V._container.assign((Real *)v_buffer.data, (Real *)v_buffer.data + _field._V.size());

        // Output u,v,p
        if (t >= output_counter * _output_freq) {
            output_vtk(timestep);
            output_counter++;
        }
        t += dt;
        timestep++;
    }
    scratch_buffer.destroy();
    rs_buffer.destroy();
    p_buffer.destroy();
    residual_buffer.destroy();
    f_buffer.destroy();
    g_buffer.destroy();
    u_buffer.destroy();
    v_buffer.destroy();
    t_buffer.destroy();
    cell_buffer.destroy();
    neighborhood_buffer.destroy();
    // Print Summary
    logger.finish();
    // Output u,v,p
    output_vtk(timestep);
    simulation.cleanup({discretization_pipeline, rs_pipeline, vel_pipeline, p_pipeline_red, p_pipeline_black,
                        residual_pipeline, p_boundary_pipeline});
}

void Case::output_vtk(int timestep, int my_rank) {
    // Create a new structured grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();

    // Create grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    Real dx = _grid.dx();
    Real dy = _grid.dy();

    Real x = _grid.domain().imin * dx;
    Real y = _grid.domain().jmin * dy;

    { y += dy; }
    { x += dx; }

    Real z = 0;
    for (int col = 0; col < _grid.domain().size_y + 1; col++) {
        x = _grid.domain().imin * dx;
        { x += dx; }
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

    // Velocity Array
    VTK_Array *Velocity = VTK_Array::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);

    // Print pressure and place ghost cells
    for (int j = 1; j < _grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _grid.domain().size_x + 1; i++) {
            Real pressure = _field.p(i, j);
            Pressure->InsertNextTuple(&pressure);

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
            vel[0] = (_field.u(i, j) + _field.u(i, j + 1)) * 0.5;
            vel[1] = (_field.v(i, j) + _field.v(i + 1, j)) * 0.5;
            Velocity->InsertNextTuple(vel);
        }
    }

    // Add Pressure to Structured Grid
    structuredGrid->GetCellData()->AddArray(Pressure);

    // Add Velocity to Structured Grid
    structuredGrid->GetPointData()->AddArray(Velocity);

    // Add Temperature to Structured Grid
    if (_calc_temp) {
        VTK_Array *Temperature = VTK_Array::New();
        Temperature->SetName("temperature");
        Temperature->SetNumberOfComponents(1);

        for (int j = 1; j < _grid.domain().size_y + 1; j++) {
            for (int i = 1; i < _grid.domain().size_x + 1; i++) {

                Real temperature = _field.t(i, j);
                Temperature->InsertNextTuple(&temperature);
            }
        }
        structuredGrid->GetCellData()->AddArray(Temperature);
    }

    // Write Grid
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();

    // Create Filename
    std::string outputname =
        _dict_name + '/' + _case_name + "_" + std::to_string(my_rank) + "." + std::to_string(timestep) + ".vtk";

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
