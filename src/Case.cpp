#include "Case.hpp"
#include "Enums.hpp"
#include "GPUSimulation.hpp"
#include <algorithm>
#include <cmath>
#include <ctime>
#include <execution>
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

#define ENABLE_GS_GPU 0

#define ENABLE_GS_CPU 0

#define ENABLE_CG_CPU 0

#define ENABLE_PRECOND 1
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

    // Check if this case uses energy equation
    if (TI != REAL_MAX) {
        _calc_temp = true;
    }

    // Prandtl number = nu / alpha
    if (pr != REAL_MAX) {
        alpha = nu / pr;
    } else if (alpha == REAL_MAX) {
        std::cerr << "Prandtl number, alpha or beta are not set, defaulting to 0\n";
        alpha = 0.0;
        beta = 0.0;
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
    _field.calc_temp = _calc_temp;
    _discretization = Discretization(domain.dx, domain.dy, gamma);
    _pressure_solver_sor = std::make_unique<SOR>(omg);
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
    auto grid_x = _grid.imaxb();
    auto grid_y = _grid.jmaxb();
    auto grid_size = grid_x * grid_y;

    GPUSimulation simulation;
    UBOData ubo_data;
    Buffer ubo_buffer;
    Buffer cell_buffer;
    Buffer neighborhood_buffer;
    Buffer u_buffer;
    Buffer v_buffer;
    Buffer f_buffer;
    Buffer g_buffer;
    Buffer t_new_buffer;
    Buffer t_old_buffer;
    Buffer rs_buffer;
    Buffer p_buffer;
    Buffer residual_buffer;
    Buffer scratch_buffer;

    Buffer a_data_buffer;
    Buffer a_offset_buffer;
    Buffer m_data_buffer;
    Buffer m_offset_buffer;

    Buffer d_buffer;
    Buffer spmv_result_buffer;
    Buffer counter_buffer;
    Buffer r_buffer;
    Buffer z_buffer;
    Buffer deltas_buffer;
    Buffer dt_buffer;
    Buffer u_boundary_matrix_buffer;
    Buffer v_boundary_matrix_buffer;
    Buffer u_rhs_buffer;
    Buffer v_rhs_buffer;
    Buffer u_row_start;
    Buffer v_row_start;
    Buffer u_col_index;
    Buffer v_col_index;
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
            BoundaryData data;
            uint32_t type = boundary->get_type();
            data.neighborhood |= type << 8;
            // data.idx = _grid.imaxb() * j + i;
            if (cell->is_border(border_position::RIGHT)) {
                data.neighborhood |= data.neighborhood | 1;
            }
            if (cell->is_border(border_position::LEFT)) {
                data.neighborhood |= data.neighborhood | 2;
            }
            if (cell->is_border(border_position::TOP)) {
                data.neighborhood |= data.neighborhood | 4;
            }
            if (cell->is_border(border_position::BOTTOM)) {
                data.neighborhood |= data.neighborhood | 8;
            }
            boundaries[j * _grid.imaxb() + i] = data;
        }
    }
    _pressure_solver_pcg =
        std::make_unique<PCG>(_grid.imaxb(), _grid.jmaxb(), _grid.dx(), _grid.dy(), _field, _grid, _boundaries);
    DiagonalSparseMatrix<Real> A_matrix_diag =
        create_diagonal_matrix(static_cast<PCG *>(_pressure_solver_pcg.get())->A, _grid.imaxb(), _grid.jmaxb(),
                               {-_grid.imaxb(), -1, 0, 1, _grid.imaxb()});
    DiagonalSparseMatrix<Real> A_precond_diag;

    UBOData data = {_grid.imaxb(),
                    _grid.jmaxb(),
                    _grid.imaxb() * _grid.jmaxb(),
                    _field._nu,
                    _field._alpha,
                    _field._beta,
                    _field._gx,
                    _field._gy,
                    _grid.dx(),
                    _grid.dy(),
                    _grid.dx() * _grid.dx(),
                    _grid.dy() * _grid.dy(),
                    1 / _grid.dx(),
                    1 / _grid.dy(),
                    _discretization._gamma,
                    _field._PI,
                    _field._UI,
                    _field._VI,
                    _field._tau};
    if (ENABLE_PRECOND) {
        A_precond_diag =
            create_preconditioner_spai(static_cast<PCG *>(_pressure_solver_pcg.get())->A, _grid.imaxb(), _grid.jmaxb());
        data.num_diags = A_precond_diag.num_diags;
    }
    simulation.init(data);
    std::vector<VkDescriptorSetLayoutBinding> set_layout_bindings;
    for (int i = 0; i < 32; i++) {
        set_layout_bindings.push_back(
            descriptor_set_layout_binding(VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, VK_SHADER_STAGE_COMPUTE_BIT, i));
    }
    simulation.create_descriptor_set_layout(set_layout_bindings);
    simulation.create_command_pool();
    simulation.create_fences();
    simulation.create_query_pool(32, VK_QUERY_TYPE_TIMESTAMP);

    Pipeline fg_pipeline = simulation.create_compute_pipeline("src/shaders/calc_fg.comp.spv", _calc_temp);
    Pipeline rs_pipeline = simulation.create_compute_pipeline("src/shaders/calc_rs.comp.spv");
    Pipeline vel_pipeline = simulation.create_compute_pipeline("src/shaders/calc_vel.comp.spv");
    Pipeline p_pipeline_red = simulation.create_compute_pipeline("src/shaders/calc_p_redblack_gs.comp.spv", 0);
    Pipeline p_pipeline_black = simulation.create_compute_pipeline("src/shaders/calc_p_redblack_gs.comp.spv", 1);
    Pipeline residual_pipeline = simulation.create_compute_pipeline("src/shaders/calc_r.comp.spv");
    Pipeline p_boundary_pipeline = simulation.create_compute_pipeline("src/shaders/boundary_p.comp.spv");
    Pipeline v_boundary_pipeline = simulation.create_compute_pipeline("src/shaders/boundary_v.comp.spv");
    Pipeline fg_boundary_pipeline = simulation.create_compute_pipeline("src/shaders/boundary_fg.comp.spv");
    Pipeline spmv_a_pipeline = simulation.create_compute_pipeline("src/shaders/spmv.comp.spv", 0);
    Pipeline saxpy_0_pipeline = simulation.create_compute_pipeline("src/shaders/saxpy.comp.spv", 0);
    Pipeline saxpy_1_pipeline = simulation.create_compute_pipeline("src/shaders/saxpy.comp.spv", 1);
    Pipeline saxpy_2_pipeline = simulation.create_compute_pipeline("src/shaders/saxpy.comp.spv", 2);
    Pipeline vec_dot_vec_0_pipeline = simulation.create_compute_pipeline("src/shaders/vec_dot_vec.comp.spv", 0);
    Pipeline vec_dot_vec_1_pipeline = simulation.create_compute_pipeline("src/shaders/vec_dot_vec.comp.spv", 1);
    Pipeline vec_dot_vec_2_pipeline;
    Pipeline spmv_m_pipeline;
    Pipeline saxpy_3_pipeline;
    Pipeline calc_t_pipeline;
    Pipeline boundary_t_branchless;
#if ENABLE_PRECOND
    vec_dot_vec_2_pipeline = simulation.create_compute_pipeline("src/shaders/vec_dot_vec.comp.spv", 2);
    spmv_m_pipeline = simulation.create_compute_pipeline("src/shaders/spmv.comp.spv", 1);
    saxpy_3_pipeline = simulation.create_compute_pipeline("src/shaders/saxpy.comp.spv", 3);
#endif
    Pipeline div_pipeline = simulation.create_compute_pipeline("src/shaders/div.comp.spv", 0);
    Pipeline div_store_pipeline = simulation.create_compute_pipeline("src/shaders/div.comp.spv", 1);
    Pipeline reduce_pipeline = simulation.create_compute_pipeline("src/shaders/reduce.comp.spv", 0);
    Pipeline inc_pipeline = simulation.create_compute_pipeline("src/shaders/increment.comp.spv");
    Pipeline negate_pipeline = simulation.create_compute_pipeline("src/shaders/negate.comp.spv");
    Pipeline min_max_uv_pipeline = simulation.create_compute_pipeline("src/shaders/min_max_uv.comp.spv");
    Pipeline reduce_u_pipeline = simulation.create_compute_pipeline("src/shaders/reduce_uv.comp.spv", 0);
    Pipeline reduce_v_pipeline = simulation.create_compute_pipeline("src/shaders/reduce_uv.comp.spv", 1);
    Pipeline calc_dt_pipeline = simulation.create_compute_pipeline("src/shaders/calc_dt.comp.spv", _calc_temp);
    Pipeline boundary_uv_branchless_pipeline =
        simulation.create_compute_pipeline("src/shaders/boundary_uv_branchless.comp.spv");
    if (_calc_temp) {
        calc_t_pipeline = simulation.create_compute_pipeline("src/shaders/calc_temp.comp.spv");
        boundary_t_branchless = simulation.create_compute_pipeline("src/shaders/boundary_t_branchless.comp.spv");
    }

    cell_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                       VK_SHARING_MODE_EXCLUSIVE, sizeof(Real) * _field.f_matrix().size(), is_fluid.data(), true);

    neighborhood_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                               VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE,
                               boundaries.size() * sizeof(BoundaryData), boundaries.data(), true);
    u_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
                    VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE,
                    _field.u_matrix().size() * sizeof(Real), _field._U._container.data(), true);
    v_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
                    VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE,
                    _field.v_matrix().size() * sizeof(Real), _field._V._container.data(), true);
    f_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                    VK_SHARING_MODE_EXCLUSIVE, _field.f_matrix().size() * sizeof(Real), _field._F._container.data(),
                    true);
    g_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                    VK_SHARING_MODE_EXCLUSIVE, _field.g_matrix().size() * sizeof(Real), _field._G._container.data(),
                    true);

    // residual_buffer.create(&simulation.context,
    //                       VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT |
    //                           VK_BUFFER_USAGE_TRANSFER_DST_BIT,
    //                       VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
    //                       VK_SHARING_MODE_EXCLUSIVE, 32 * sizeof(Real));
    residual_buffer.create(
        &simulation.context,
        VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE, ceil(grid_size / 32) * sizeof(Real));

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

    if (ENABLE_PRECOND) {
        m_data_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                             VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE,
                             A_precond_diag.data.size() * sizeof(Real), A_precond_diag.data.data(), true);
        m_offset_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                               VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE,
                               A_precond_diag.offsets.size() * sizeof(int), A_precond_diag.offsets.data(), true);
        z_buffer.create(
            &simulation.context,
            VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE, _field.p_matrix().size() * sizeof(Real));
    }

    a_data_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                         VK_SHARING_MODE_EXCLUSIVE, A_matrix_diag.data.size() * sizeof(Real), A_matrix_diag.data.data(),
                         true);
    a_offset_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                           VK_SHARING_MODE_EXCLUSIVE, 5 * sizeof(int), A_matrix_diag.offsets.data(), true);

    d_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                    VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE,
                    _field.p_matrix().size() * sizeof(Real));
    spmv_result_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                              VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE,
                              _field.p_matrix().size() * sizeof(Real));
    counter_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                          VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE, sizeof(int));
    /*counter_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                         VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                         VK_SHARING_MODE_EXCLUSIVE, sizeof(int));
    *(int *)counter_buffer.data = 0;*/
    deltas_buffer.create(&simulation.context,
                         VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT |
                             VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                         VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE, 2 * sizeof(Real));
    r_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                    VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE,
                    _field.p_matrix().size() * sizeof(Real));

    /*  dt_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                       VK_SHARING_MODE_EXCLUSIVE, sizeof(Real));*/
    dt_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                     VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                     VK_SHARING_MODE_EXCLUSIVE, sizeof(Real));
    Buffer t_rhs_buffer;
    Buffer t_row_start;
    Buffer t_col_index;
    Buffer t_boundary_matrix_buffer;

    auto pcg_solver = (PCG *)_pressure_solver_pcg.get();
    auto u_matrix_data = pcg_solver->U_fixed.value.data();
    auto u_matrix_size = pcg_solver->U_fixed.value.size();
    auto v_matrix_data = pcg_solver->V_fixed.value.data();
    auto v_matrix_size = pcg_solver->V_fixed.value.size();
    auto u_row_start_data = pcg_solver->U_fixed.rowstart.data();
    auto u_row_start_size = pcg_solver->U_fixed.rowstart.size();
    auto u_col_idx_data = pcg_solver->U_fixed.colindex.data();
    auto u_col_idx_size = pcg_solver->U_fixed.colindex.size();
    auto v_row_start_data = pcg_solver->V_fixed.rowstart.data();
    auto v_row_start_size = pcg_solver->V_fixed.rowstart.size();
    auto v_col_idx_data = pcg_solver->V_fixed.colindex.data();
    auto v_col_idx_size = pcg_solver->V_fixed.colindex.size();
    auto u_rhs_data = pcg_solver->U_RHS.data();
    auto u_rhs_size = pcg_solver->U_RHS.size();
    auto v_rhs_data = pcg_solver->V_RHS.data();
    auto v_rhs_size = pcg_solver->V_RHS.size();

    u_boundary_matrix_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                                    VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE,
                                    u_matrix_size * sizeof(Real), u_matrix_data, true);
    v_boundary_matrix_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                                    VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE,
                                    v_matrix_size * sizeof(Real), v_matrix_data, true);
    u_rhs_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                        VK_SHARING_MODE_EXCLUSIVE, u_rhs_size * sizeof(Real), u_rhs_data, true);
    v_rhs_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                        VK_SHARING_MODE_EXCLUSIVE, v_rhs_size * sizeof(int), v_rhs_data, true);
    u_row_start.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                       VK_SHARING_MODE_EXCLUSIVE, u_row_start_size * sizeof(int), u_row_start_data, true);
    v_row_start.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                       VK_SHARING_MODE_EXCLUSIVE, v_row_start_size * sizeof(int), v_row_start_data, true);
    u_col_index.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                       VK_SHARING_MODE_EXCLUSIVE, u_col_idx_size * sizeof(int), u_col_idx_data, true);
    v_col_index.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                       VK_SHARING_MODE_EXCLUSIVE, v_col_idx_size * sizeof(int), v_col_idx_data, true);

    if (_calc_temp) {
        auto t_matrix_data = pcg_solver->T_fixed.value.data();
        auto t_matrix_size = pcg_solver->T_fixed.value.size();
        auto t_row_start_data = pcg_solver->T_fixed.rowstart.data();
        auto t_row_start_size = pcg_solver->T_fixed.rowstart.size();
        auto t_col_idx_data = pcg_solver->T_fixed.colindex.data();
        auto t_col_idx_size = pcg_solver->T_fixed.colindex.size();
        auto t_rhs_data = pcg_solver->T_RHS.data();
        auto t_rhs_size = pcg_solver->T_RHS.size();
        t_new_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
                            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE,
                            _field.t_matrix().size() * sizeof(Real), _field._T._container.data(), true);
        t_old_buffer.create(
            &simulation.context,
            VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE, _field.t_matrix().size() * sizeof(Real));
        t_boundary_matrix_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                                        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE,
                                        t_matrix_size * sizeof(Real), t_matrix_data, true);
        t_rhs_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE, t_rhs_size * sizeof(Real),
                            t_rhs_data, true);
        t_row_start.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                           VK_SHARING_MODE_EXCLUSIVE, t_row_start_size * sizeof(int), t_row_start_data, true);
        t_col_index.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                           VK_SHARING_MODE_EXCLUSIVE, t_col_idx_size * sizeof(int), t_col_idx_data, true);
    }

    simulation.push_descriptors({
        {u_buffer, 0}, {v_buffer, 1}, {f_buffer, 2}, {g_buffer, 3}, {cell_buffer, 5}, {rs_buffer, 6}, {p_buffer, 7},
            {residual_buffer, 8}, {neighborhood_buffer, 9}, {a_data_buffer, 10}, {a_offset_buffer, 11}, {d_buffer, 12},
            {spmv_result_buffer, 13}, {residual_buffer, 14}, {r_buffer, 15}, {p_buffer, 16},
#if ENABLE_PRECOND
            {z_buffer, 17}, {m_data_buffer, 18}, {m_offset_buffer, 19},
#endif
            {dt_buffer, 21}, {deltas_buffer, 30}, {counter_buffer, 31}, {u_boundary_matrix_buffer, 0, 1},
            {v_boundary_matrix_buffer, 1, 1}, {u_rhs_buffer, 2, 1}, {v_rhs_buffer, 3, 1}, {u_row_start, 4, 1},
            {v_row_start, 5, 1}, {u_col_index, 6, 1}, {v_col_index, 7, 1},
    });
    if (_calc_temp) {
        simulation.push_descriptors({{t_old_buffer, 22},
                                     {t_new_buffer, 23},
                                     {t_boundary_matrix_buffer, 8, 1},
                                     {t_rhs_buffer, 9, 1},
                                     {t_row_start, 10, 1},
                                     {t_col_index, 11, 1}});
    }

    auto record_simulation_step = [&](int command_idx) {
        simulation.begin_recording(command_idx, 0);
        vkCmdResetQueryPool(simulation.context.command_buffer[command_idx], simulation.query_pool, 0, 6);
        simulation.write_timestamp(command_idx, 0);
        simulation.record_command_buffer(min_max_uv_pipeline, command_idx, 1024, 1, grid_size, 1);
        uv_max(simulation, residual_buffer, counter_buffer, min_max_uv_pipeline, reduce_u_pipeline, reduce_v_pipeline,
               command_idx, grid_size);
        barrier(simulation, residual_buffer, command_idx);
        simulation.record_command_buffer(calc_dt_pipeline, command_idx, 1, 1, 1, 1);
        barrier(simulation, dt_buffer, command_idx);
        simulation.record_command_buffer(fg_pipeline, command_idx, 32, 32, grid_x, grid_y);
        std::array<VkBufferMemoryBarrier, 2> fg_barriers = {
            buffer_barrier(f_buffer.handle, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_WRITE_BIT),
            buffer_barrier(g_buffer.handle, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_WRITE_BIT)

        };
        vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                             VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 2, fg_barriers.data(), 0, 0);
        if (_calc_temp) {
            t_new_buffer.copy(t_old_buffer, command_idx, false);
            simulation.record_command_buffer(boundary_t_branchless, command_idx, 1024, 1, grid_size, 1);
            barrier(simulation, t_new_buffer, command_idx);
            simulation.record_command_buffer(calc_t_pipeline, command_idx, 32, 32, grid_x, grid_y);
            barrier(simulation, t_new_buffer, command_idx);
        }
        // Compute F & G and enforce boundary conditions
        simulation.record_command_buffer(fg_boundary_pipeline, command_idx, 32, 32, grid_x, grid_y);

        std::array<VkBufferMemoryBarrier, 2> fg_barriers_read = {
            buffer_barrier(f_buffer.handle, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT),
            buffer_barrier(g_buffer.handle, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT)};
        vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                             VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 2, fg_barriers_read.data(), 0, 0);
        simulation.record_command_buffer(rs_pipeline, command_idx, 32, 32, grid_x, grid_y);
        vkCmdFillBuffer(simulation.context.command_buffer[command_idx], p_buffer.handle, 0,
                        _grid.imaxb() * _grid.jmaxb() * sizeof(Real), 0);
        std::vector<Real> solution(_field.p_matrix().size(), 0);
        VkBufferMemoryBarrier fill_barrier = buffer_barrier(p_buffer.handle, VK_ACCESS_TRANSFER_WRITE_BIT,
                                                            VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_SHADER_READ_BIT);
        vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_TRANSFER_BIT,
                             VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &fill_barrier, 0, 0);
        barrier(simulation, rs_buffer, command_idx);

        rs_buffer.copy(r_buffer, command_idx, false);
        barrier(simulation, r_buffer, command_idx);

        if (!ENABLE_PRECOND) {
            // z_buffer.copy(d_buffer, command_idx, false);
            rs_buffer.copy(d_buffer, command_idx, false);

            vec_dp(simulation, residual_buffer, counter_buffer, vec_dot_vec_1_pipeline, reduce_pipeline, command_idx,
                   grid_size);
        } else {
            simulation.record_command_buffer(spmv_m_pipeline, command_idx, 1024, 1, _grid.imaxb() * _grid.jmaxb(), 1);
            barrier(simulation, z_buffer, command_idx);
            z_buffer.copy(d_buffer, command_idx, false);
            vec_dp(simulation, residual_buffer, counter_buffer, vec_dot_vec_2_pipeline, reduce_pipeline, command_idx,
                   grid_size);
        }
        barrier(simulation, residual_buffer, command_idx);
        residual_buffer.copy(scratch_buffer, command_idx, false, 0, 0, sizeof(Real));
        residual_buffer.copy(deltas_buffer, command_idx, false, 0, 0, sizeof(Real));
        // Calculate delta_new(r_norm) = r_t dot r
        simulation.write_timestamp(command_idx, 1);
        simulation.end_recording(command_idx);
    };

    auto record_conjugate_gradient_solver = [&](int command_idx = 0) {
        simulation.begin_recording(command_idx, VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT);
        vkCmdResetQueryPool(simulation.context.command_buffer[command_idx], simulation.query_pool, 0, 6);
        simulation.write_timestamp(command_idx, 2);
        // q <- A *d
        simulation.record_command_buffer(spmv_a_pipeline, command_idx, 1024, 1, _grid.imaxb() * _grid.jmaxb(), 1);
        barrier(simulation, spmv_result_buffer, command_idx, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT);
        // Store d^T dot q scalar in the residual buffer
        vec_dp(simulation, residual_buffer, counter_buffer, vec_dot_vec_0_pipeline, reduce_pipeline, command_idx,
               grid_size);
        barrier(simulation, residual_buffer, command_idx);

        // alpha = delta_new / dt_dotq
        scalar_div(simulation, div_pipeline, command_idx);
        barrier(simulation, residual_buffer, command_idx);

        // x <- x + alpha * d
        vec_saxpy(simulation, saxpy_0_pipeline, command_idx, grid_size);
        //barrier(simulation, p_buffer, command_idx);
        // r <- r - alpha * q
        vec_saxpy(simulation, saxpy_1_pipeline, command_idx, grid_size);
        barrier(simulation, r_buffer, command_idx, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT);
       
        // simulation.record_command_buffer(inc_pipeline, command_idx, 1, 1, 1, 1);
        // barrier(simulation, counter_buffer, command_idx);
        if (ENABLE_PRECOND) {
            // Preconditioning
            // z <- M *r
            simulation.record_command_buffer(spmv_m_pipeline, command_idx, 1024, 1, _grid.imaxb() * _grid.jmaxb(), 1);
            barrier(simulation, z_buffer, command_idx, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT);

            vec_dp(simulation, residual_buffer, counter_buffer, vec_dot_vec_2_pipeline, reduce_pipeline, command_idx,
                   grid_size);
        } else {
            vec_dp(simulation, residual_buffer, counter_buffer, vec_dot_vec_1_pipeline, reduce_pipeline, command_idx,
                   grid_size);
        }

        barrier(simulation, residual_buffer, command_idx);
        /*  delta_old = delta_new;
          delta_new = vec_dp(simulation, r_buffer, r_buffer, residual_buffer, scratch_buffer, vec_dot_vec_pipeline,
                             reduce_pipeline);
          Real beta = delta_new / delta_old;*/
        scalar_div(simulation, div_store_pipeline, command_idx);
        barrier(simulation, residual_buffer, command_idx, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT);
        if (ENABLE_PRECOND) {
            // d <- z + beta * d;
            vec_saxpy(simulation, saxpy_3_pipeline, command_idx, grid_size);
        } else {
            // d <- r + beta * d;
            vec_saxpy(simulation, saxpy_2_pipeline, command_idx, grid_size);
        }

        // simulation.record_command_buffer(inc_pipeline, command_idx, 1, 1, 1, 1);
        // barrier(simulation, counter_buffer, command_idx);
       
      
        simulation.write_timestamp(command_idx, 3);
        simulation.end_recording(command_idx);
    };

    auto record_post_pressure = [&](int command_idx) {
        simulation.begin_recording(command_idx, 0);
        vkCmdResetQueryPool(simulation.context.command_buffer[command_idx], simulation.query_pool, 0, 6);
        simulation.write_timestamp(command_idx, 4);
        simulation.record_command_buffer(negate_pipeline, command_idx, 1024, 1, grid_size, 1);
        barrier(simulation, p_buffer, command_idx);
        simulation.record_command_buffer(vel_pipeline, command_idx, 32, 32, grid_x, grid_y);
        std::array<VkBufferMemoryBarrier, 2> uv_barriers = {
            buffer_barrier(u_buffer.handle, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_WRITE_BIT),
            buffer_barrier(v_buffer.handle, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_WRITE_BIT)};
        vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                             VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 2, uv_barriers.data(), 0, 0);
        // simulation.record_command_buffer(v_boundary_pipeline);
        simulation.record_command_buffer(boundary_uv_branchless_pipeline, command_idx, 1024, 1, grid_size, 1);
        simulation.write_timestamp(command_idx, 5);
        simulation.end_recording(command_idx);
    };
    record_simulation_step(0);
    record_conjugate_gradient_solver(1);
    record_post_pressure(2);
    while (t < _t_end) {
        auto c_start = std::chrono::high_resolution_clock::now();

        // std::cout << t << "/" << _t_end;
        // Select dt
        // dt = _field.calculate_dt(_grid, _calc_temp);

        // Enforce velocity boundary conditions
        /* for (auto &boundary : _boundaries) {
             boundary->enforce_uv(_field);
         } */

        // if (_calc_temp) {
        //    // Enforce temperature boundary conditions
        //    for (const auto &boundary : _boundaries) {
        //        boundary->enforce_t(_field);
        //    }
        //    // Compute temperatures
        //    _field.calculate_temperatures(_grid);
        //}

        // Compute F & G and enforce boundary conditions
        //_field.calculate_fluxes(_grid, _calc_temp);

        // for (const auto &boundary : _boundaries) {
        //    boundary->enforce_fg(_field);
        //}
        //// Set RHS of PPE
        ////_field.calculate_rs(_grid);

        // Perform pressure solve
        uint32_t it = 0;
        Real res = REAL_MAX;
#if ENABLE_CG_CPU || ENABLE_GS_CPU
        rs_buffer.copy(scratch_buffer);
        _field._RS._container.assign((Real *)scratch_buffer.data, (Real *)scratch_buffer.data + _field._RS.size());
#endif
#if ENABLE_CG_CPU

        ////// CG-CPU
        res = _pressure_solver_pcg->solve(_field, _grid, _boundaries, _max_iter, _tolerance);
#endif

#if ENABLE_GS_CPU
        // GS-CPU
        while (it < _max_iter && res > _tolerance) {
            res = _pressure_solver_sor->solve(_field, _grid, _boundaries, _max_iter, _tolerance);
            // Enforce boundary conditions
            for (const auto &boundary : _boundaries) {
                boundary->enforce_p(_field);
            }
            it++;
        }
#endif

#if ENABLE_GS_GPU
        it = 0;

        // Red-Black Gauss Seidel-GPU
        // std::cout << " ------- " << res;
        while (it < _max_iter && res > _tolerance) {
            simulation.begin_recording();
            vkCmdFillBuffer(simulation.context.command_buffer, residual_buffer.handle, 0, 32 * sizeof(Real), 0);
            simulation.record_command_buffer(p_pipeline_red);
            auto p_barrier = buffer_barrier(p_buffer.handle, VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_SHADER_READ_BIT,
                                            VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_SHADER_READ_BIT);
            vkCmdPipelineBarrier(simulation.context.command_buffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                                 VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &p_barrier, 0, 0);
            simulation.record_command_buffer(p_pipeline_black);
            vkCmdPipelineBarrier(simulation.context.command_buffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                                 VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &p_barrier, 0, 0);
            VkBufferMemoryBarrier fill_barrier =
                buffer_barrier(residual_buffer.handle, VK_ACCESS_TRANSFER_WRITE_BIT, VK_ACCESS_SHADER_WRITE_BIT);
            vkCmdPipelineBarrier(simulation.context.command_buffer, VK_PIPELINE_STAGE_TRANSFER_BIT,
                                 VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &fill_barrier, 0, 0);
            simulation.record_command_buffer(residual_pipeline);
            simulation.record_command_buffer(p_boundary_pipeline);
            simulation.end_recording();
            simulation.run_command_buffer();
            residual_buffer.copy(scratch_buffer);
            /*  std::vector<Real> test(32);
              test.assign((Real *)scratch_buffer.data, (Real *)scratch_buffer.data + 32);
              Real test2 = test[0];*/
            res = std::reduce(std::execution::par, (Real *)scratch_buffer.data, (Real *)scratch_buffer.data + 32);
            res = std::sqrt(res / _grid.fluid_cells().size());
            it++;
        }
#endif

        simulation.run_command_buffer(0);
        uint64_t timestamps[6] = {};
        simulation.get_query_results(COUNTOF(timestamps), timestamps);
        double begin = double(timestamps[0]) * simulation.props.limits.timestampPeriod * 1e-6;
        double end = double(timestamps[1]) * simulation.props.limits.timestampPeriod * 1e-6;

        /*  auto dtgpu = *(Real *)dt_buffer.data;
          std::vector<Real> umaxbuf(32);
          std::vector<Real> vmaxbuf(32);

          umaxbuf.assign((Real *)u_max_buffer.data, (Real *)u_max_buffer.data + 32);
          vmaxbuf.assign((Real *)v_max_buffer.data, (Real *)v_max_buffer.data + 32);
           dt_buffer.upload(sizeof(float), 0, &dt, &scratch_buffer);*/
        Real delta_new = *(Real *)scratch_buffer.data;
        Real delta_old = delta_new;
        Real delta_zero = delta_new;
        Real cond = _tolerance * _tolerance * delta_zero;
        double pressure_time = 0;
        std::chrono::nanoseconds ptimecpu(0);
        uint8_t sem_idx = 0;
        VkBufferCopy copy_region;
        copy_region.srcOffset = 0;
        copy_region.dstOffset = 0;
        copy_region.size = deltas_buffer.size;
        while (it < _max_iter && delta_new > cond) {
            for (int i = 0; i < 35; i++) {
                VkSubmitInfo submit_info = {};
                submit_info.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
                submit_info.commandBufferCount = 1;
                submit_info.pCommandBuffers = &simulation.context.command_buffer[1];
                submit_info.pWaitSemaphores = i == 0 ? 0 : &simulation.semaphores[sem_idx]; 
                submit_info.pSignalSemaphores = &simulation.semaphores[sem_idx ^ 1]; 
                VK_CHECK(vkQueueSubmit(simulation.context.compute_queue, 1, &submit_info, simulation.fences[i]));
                sem_idx ^= 1;
            }
            it += 35;
            vkWaitForFences(simulation.context.device, simulation.fences.size(), simulation.fences.data(), true,
                            100000000000);
            vkResetFences(simulation.context.device, simulation.fences.size(), simulation.fences.data());
            simulation.get_query_results(COUNTOF(timestamps), timestamps);
            double pbegin = double(timestamps[2]) * simulation.props.limits.timestampPeriod * 1e-6;
            double pend = double(timestamps[3]) * simulation.props.limits.timestampPeriod * 1e-6;
            pressure_time += pend - pbegin;
            simulation.begin_recording(3, 0);
            barrier(simulation, deltas_buffer, 3, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT);
            vkCmdCopyBuffer(simulation.context.command_buffer[3], deltas_buffer.handle, scratch_buffer.handle,
                            1, &copy_region);
            simulation.end_recording(3);
            simulation.run_command_buffer(3);
            delta_new = *(Real *)scratch_buffer.data;
        }
        if (it == _max_iter) {
            // std::cout << " ------ " << it << " " << res;
            logger.max_iter_warning();
        }
        // Output current timestep information
        logger.write_log(timestep, t, it, _max_iter, res);

        // Compute u^(n+1) & v^(n+1)
        //_field.calculate_velocities(_grid);

        simulation.run_command_buffer(2);
        simulation.get_query_results(COUNTOF(timestamps), timestamps);
        double post_begin = double(timestamps[4]) * simulation.props.limits.timestampPeriod * 1e-6;
        double post_end = double(timestamps[5]) * simulation.props.limits.timestampPeriod * 1e-6;

        u_buffer.copy(scratch_buffer, 3);
        _field._U._container.assign((Real *)scratch_buffer.data, (Real *)scratch_buffer.data + _field._U.size());
        v_buffer.copy(scratch_buffer, 3);
        _field._V._container.assign((Real *)scratch_buffer.data, (Real *)scratch_buffer.data + _field._V.size());
        p_buffer.copy(scratch_buffer, 3);
        _field._P._container.assign((Real *)scratch_buffer.data, (Real *)scratch_buffer.data + _field._P.size());
        if (_calc_temp) {
            t_new_buffer.copy(scratch_buffer, 3);
            _field._T._container.assign((Real *)scratch_buffer.data, (Real *)scratch_buffer.data + _field._T.size());
        }

        // Output u,v,p
        if (t >= output_counter * _output_freq) {
            output_vtk(timestep);
            output_counter++;
        }
        dt = *(Real *)dt_buffer.data;
        t += dt;
        timestep++;
        auto c_end = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start);
        // auto time = std::chrono::duration_cast<std::chrono::milliseconds>(ptimecpu);
        printf("\rIter: %d, first pass time %.2f ms, pressure time %.2f ms, post pressure %.2f ms, CPU total: %ld ms", it,
               end - begin, pressure_time, post_end - post_begin, time.count());
        // Print progress bar
        logger.progress_bar(t, _t_end);
    }
    scratch_buffer.destroy();
    rs_buffer.destroy();
    p_buffer.destroy();
    residual_buffer.destroy();
    f_buffer.destroy();
    g_buffer.destroy();
    u_buffer.destroy();
    v_buffer.destroy();
    ubo_buffer.destroy();
    cell_buffer.destroy();
    a_data_buffer.destroy();
    a_offset_buffer.destroy();
#if ENABLE_PRECOND
    m_data_buffer.destroy();
    m_offset_buffer.destroy();
    z_buffer.destroy();
#endif
    d_buffer.destroy();
    spmv_result_buffer.destroy();
    counter_buffer.destroy();
    r_buffer.destroy();
    deltas_buffer.destroy();
    neighborhood_buffer.destroy();
    dt_buffer.destroy();
    u_boundary_matrix_buffer.destroy();
    v_boundary_matrix_buffer.destroy();
    u_rhs_buffer.destroy();
    v_rhs_buffer.destroy();
    u_row_start.destroy();
    v_row_start.destroy();
    u_col_index.destroy();
    v_col_index.destroy();
    if (_calc_temp) {
        t_new_buffer.destroy();
        t_old_buffer.destroy();
        t_rhs_buffer.destroy();
        t_row_start.destroy();
        t_col_index.destroy();
        t_boundary_matrix_buffer.destroy();
    }
    // Print Summary
    logger.finish();
    // Output u,v,p
    output_vtk(timestep);
    std::vector<Pipeline> pipelines_to_destroy = {
        fg_pipeline,
        rs_pipeline,
        vel_pipeline,
        p_pipeline_red,
        p_pipeline_black,
        residual_pipeline,
        p_boundary_pipeline,
        v_boundary_pipeline,
        fg_boundary_pipeline,
        spmv_a_pipeline,
        saxpy_0_pipeline,
        saxpy_1_pipeline,
        saxpy_2_pipeline,
        reduce_pipeline,
        vec_dot_vec_0_pipeline,
        vec_dot_vec_1_pipeline,
#if ENABLE_PRECOND
        vec_dot_vec_2_pipeline,
        spmv_m_pipeline,
        saxpy_3_pipeline,
#endif
        div_pipeline,
        div_store_pipeline,
        inc_pipeline,
        negate_pipeline,
        min_max_uv_pipeline,
        reduce_u_pipeline,
        reduce_v_pipeline,
        calc_dt_pipeline,
        boundary_uv_branchless_pipeline
    };
    if (_calc_temp) {
        pipelines_to_destroy.push_back(calc_t_pipeline);
        pipelines_to_destroy.push_back(boundary_t_branchless);
    }
    simulation.cleanup(pipelines_to_destroy);
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
