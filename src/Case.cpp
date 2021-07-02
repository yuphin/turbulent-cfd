#include "Case.hpp"


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

#define ENABLE_PRECOND 0
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
    int pressure_solver = 0;
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
                if (var == "iproc") file >> params.iproc;
                if (var == "jproc") file >> params.jproc;
                if (var == "solver") file >> pressure_solver;
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
    _solver = std::make_unique<CPUSolver>();
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
        _solver->_calc_temp = true;
    }

    // Prandtl number = nu / alpha
    if (pr != REAL_MAX) {
        alpha = nu / pr;
    } else if (alpha == REAL_MAX) {
        if (params.world_rank == 0 && _solver->_calc_temp) {
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

    global_size_x = imax;
    global_size_y = jmax;
    std::vector<std::vector<int>> global_geometry;
    if (_geom_name.compare("NONE")) {

        global_geometry = parse_geometry_file(_geom_name, imax, jmax);
    } else {
        global_geometry = build_lid_driven_cavity(imax, jmax);
    }
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
   

    _solver->_grid = Grid(_geom_name, domain, local_geometry);
    _solver->_field = Fields(nu, dt, tau, domain.size_x, domain.size_y, UI, VI, PI, TI, alpha, beta, GX, GY);
    _solver->_discretization = Discretization(domain.dx, domain.dy, gamma);
    _solver->_max_iter = itermax;
    _solver->_tolerance = eps;

    // Construct boundaries
    if (!_solver->_grid.noslip_wall_cells().empty()) {
        _solver->_boundaries.push_back(
            std::make_unique<NoSlipWallBoundary>(&_solver->_grid.noslip_wall_cells(), wall_vels, wall_temps));
    }
    if (!_solver->_grid.freeslip_wall_cells().empty()) {
        _solver->_boundaries.push_back(
            std::make_unique<FreeSlipWallBoundary>(&_solver->_grid.freeslip_wall_cells(), wall_vels, wall_temps));
    }
    if (!_solver->_grid.outlet_cells().empty()) {
        _solver->_boundaries.push_back(std::make_unique<OutletBoundary>(&_solver->_grid.outlet_cells()));
    }
    if (!_solver->_grid.inlet_cells().empty()) {
        _solver->_boundaries.push_back(
            std::make_unique<InletBoundary>(&_solver->_grid.inlet_cells(), inlet_Us, inlet_Vs, inlet_Ts, DP));
    }

    switch (pressure_solver) {
    case 0: {

        _solver->_pressure_solver = std::make_unique<SOR>(omg);
        break;
    }
    case 1: {
        _solver->_pressure_solver =
            std::make_unique<PCG>(_solver->_grid.imaxb(), _solver->_grid.jmaxb(), _solver->_grid.dx(),
                                  _solver->_grid.dy(), _solver->_field, _solver->_grid, _solver->_boundaries);
        break;
    }
    default:
        break;
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
    }
    _solver->solve_post_pressure();
}

// void Simulation::simulate(Params &params) {
//    Real t = 0.0;
//    Real dt = _field.dt();
//    uint32_t timestep = 0;
//    Real output_counter = 0.0;
//    auto grid_x = _grid.imaxb();
//    auto grid_y = _grid.jmaxb();
//    auto grid_size = grid_x * grid_y;
//    while (t < _t_end) {
//        // Print progress bar
//        if (params.world_rank == 0) logger.progress_bar(t, _t_end);
//
//        // Select dt
//        dt = _field.calculate_dt(_grid, _calc_temp);
//
//      
//       
//
//      
//
//        auto record_simulation_step = [&](int command_idx) {
//       
//        };
//
//        auto record_conjugate_gradient_solver = [&](int command_idx = 0) {
//           
//        };
//
//        auto record_post_pressure = [&](int command_idx) {
//         
//        };
//        record_simulation_step(0);
//        record_conjugate_gradient_solver(1);
//        record_post_pressure(2);
//        while (t < _t_end) {
//            auto c_start = std::chrono::high_resolution_clock::now();
//
//            // std::cout << t << "/" << _t_end;
//            // Select dt
//            // dt = _field.calculate_dt(_grid, _calc_temp);
//
//            // Enforce velocity boundary conditions
//            /* for (auto &boundary : _boundaries) {
//                 boundary->enforce_uv(_field);
//             } */
//
//            // if (_calc_temp) {
//            //    // Enforce temperature boundary conditions
//            //    for (const auto &boundary : _boundaries) {
//            //        boundary->enforce_t(_field);
//            //    }
//            //    // Compute temperatures
//            //    _field.calculate_temperatures(_grid);
//            //}
//
//            // Compute F & G and enforce boundary conditions
//            //_field.calculate_fluxes(_grid, _calc_temp);
//
//            // for (const auto &boundary : _boundaries) {
//            //    boundary->enforce_fg(_field);
//            //}
//            //// Set RHS of PPE
//            ////_field.calculate_rs(_grid);
//
//            // Perform pressure solve
//            uint32_t it = 0;
//            Real res = REAL_MAX;
//#if ENABLE_CG_CPU || ENABLE_GS_CPU
//            rs_buffer.copy(scratch_buffer);
//            _field._RS._container.assign((Real *)scratch_buffer.data, (Real *)scratch_buffer.data +
//            _field._RS.size());
//#endif
//#if ENABLE_CG_CPU
//
//            ////// CG-CPU
//            res = _pressure_solver_pcg->solve(_field, _grid, _boundaries, _max_iter, _tolerance);
//#endif
//
//#if ENABLE_GS_CPU
//            // GS-CPU
//            while (it < _max_iter && res > _tolerance) {
//                res = _pressure_solver->solve(_field, _grid, _boundaries, params);
//                // Enforce boundary conditions
//                for (const auto &boundary : _boundaries) {
//                    boundary->enforce_p(_field);
//                }
//                it++;
//            }
//#endif
//
//#if ENABLE_GS_GPU
//            it = 0;
//
//            // Red-Black Gauss Seidel-GPU
//            // std::cout << " ------- " << res;
//            while (it < _max_iter && res > _tolerance) {
//                simulation.begin_recording();
//                vkCmdFillBuffer(simulation.context.command_buffer, residual_buffer.handle, 0, 32 * sizeof(Real), 0);
//                simulation.record_command_buffer(p_pipeline_red);
//                auto p_barrier = buffer_barrier(p_buffer.handle, VK_ACCESS_SHADER_WRITE_BIT |
//                VK_ACCESS_SHADER_READ_BIT,
//                                                VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_SHADER_READ_BIT);
//                vkCmdPipelineBarrier(simulation.context.command_buffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
//                                     VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &p_barrier, 0, 0);
//                simulation.record_command_buffer(p_pipeline_black);
//                vkCmdPipelineBarrier(simulation.context.command_buffer, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
//                                     VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &p_barrier, 0, 0);
//                VkBufferMemoryBarrier fill_barrier =
//                    buffer_barrier(residual_buffer.handle, VK_ACCESS_TRANSFER_WRITE_BIT, VK_ACCESS_SHADER_WRITE_BIT);
//                vkCmdPipelineBarrier(simulation.context.command_buffer, VK_PIPELINE_STAGE_TRANSFER_BIT,
//                                     VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &fill_barrier, 0, 0);
//                simulation.record_command_buffer(residual_pipeline);
//                simulation.record_command_buffer(p_boundary_pipeline);
//                simulation.end_recording();
//                simulation.run_command_buffer();
//                residual_buffer.copy(scratch_buffer);
//                /*  std::vector<Real> test(32);
//                  test.assign((Real *)scratch_buffer.data, (Real *)scratch_buffer.data + 32);
//                  Real test2 = test[0];*/
//                res = std::reduce(std::execution::par, (Real *)scratch_buffer.data, (Real *)scratch_buffer.data + 32);
//                res = std::sqrt(res / _grid.fluid_cells().size());
//                it++;
//            }
//#endif
//
//            simulation.run_command_buffer(0);
//            uint64_t timestamps[6] = {};
//            simulation.get_query_results(COUNTOF(timestamps), timestamps);
//            double begin = double(timestamps[0]) * simulation.props.limits.timestampPeriod * 1e-6;
//            double end = double(timestamps[1]) * simulation.props.limits.timestampPeriod * 1e-6;
//
//            /*  auto dtgpu = *(Real *)dt_buffer.data;
//              std::vector<Real> umaxbuf(32);
//              std::vector<Real> vmaxbuf(32);
//
//              umaxbuf.assign((Real *)u_max_buffer.data, (Real *)u_max_buffer.data + 32);
//              vmaxbuf.assign((Real *)v_max_buffer.data, (Real *)v_max_buffer.data + 32);
//               dt_buffer.upload(sizeof(float), 0, &dt, &scratch_buffer);*/
//            Real delta_new = *(Real *)scratch_buffer.data;
//            Real delta_old = delta_new;
//            Real delta_zero = delta_new;
//            Real cond = _tolerance * _tolerance * delta_zero;
//            double pressure_time = 0;
//            std::chrono::nanoseconds ptimecpu(0);
//            uint8_t sem_idx = 0;
//            VkBufferCopy copy_region;
//            copy_region.srcOffset = 0;
//            copy_region.dstOffset = 0;
//            copy_region.size = deltas_buffer.size;
//            while (it < _max_iter && delta_new > cond) {
//                for (int i = 0; i < 35; i++) {
//                    VkSubmitInfo submit_info = {};
//                    submit_info.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
//                    submit_info.commandBufferCount = 1;
//                    submit_info.pCommandBuffers = &simulation.context.command_buffer[1];
//                    submit_info.pWaitSemaphores = i == 0 ? 0 : &simulation.semaphores[sem_idx];
//                    submit_info.pSignalSemaphores = &simulation.semaphores[sem_idx ^ 1];
//                    VK_CHECK(vkQueueSubmit(simulation.context.compute_queue, 1, &submit_info, simulation.fences[i]));
//                    sem_idx ^= 1;
//                }
//                it += 35;
//                vkWaitForFences(simulation.context.device, simulation.fences.size(), simulation.fences.data(), true,
//                                100000000000);
//                vkResetFences(simulation.context.device, simulation.fences.size(), simulation.fences.data());
//                simulation.get_query_results(COUNTOF(timestamps), timestamps);
//                double pbegin = double(timestamps[2]) * simulation.props.limits.timestampPeriod * 1e-6;
//                double pend = double(timestamps[3]) * simulation.props.limits.timestampPeriod * 1e-6;
//                pressure_time += pend - pbegin;
//                simulation.begin_recording(3, 0);
//                barrier(simulation, deltas_buffer, 3, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT);
//                vkCmdCopyBuffer(simulation.context.command_buffer[3], deltas_buffer.handle, scratch_buffer.handle, 1,
//                                &copy_region);
//                simulation.end_recording(3);
//                simulation.run_command_buffer(3);
//                delta_new = *(Real *)scratch_buffer.data;
//            }
//            if (params.world_rank == 0 && it == _max_iter) {
//                // std::cout << " ------ " << it << " " << res;
//                logger.max_iter_warning();
//            }
//            // Output current timestep information
//            logger.write_log(timestep, t, dt, it, _max_iter, res);
//
//            // Compute u^(n+1) & v^(n+1)
//            //_field.calculate_velocities(_grid);
//
//            simulation.run_command_buffer(2);
//            simulation.get_query_results(COUNTOF(timestamps), timestamps);
//            double post_begin = double(timestamps[4]) * simulation.props.limits.timestampPeriod * 1e-6;
//            double post_end = double(timestamps[5]) * simulation.props.limits.timestampPeriod * 1e-6;
//
//            u_buffer.copy(scratch_buffer, 3);
//            _field._U._container.assign((Real *)scratch_buffer.data, (Real *)scratch_buffer.data + _field._U.size());
//            v_buffer.copy(scratch_buffer, 3);
//            _field._V._container.assign((Real *)scratch_buffer.data, (Real *)scratch_buffer.data + _field._V.size());
//            p_buffer.copy(scratch_buffer, 3);
//            _field._P._container.assign((Real *)scratch_buffer.data, (Real *)scratch_buffer.data + _field._P.size());
//            if (_calc_temp) {
//                t_new_buffer.copy(scratch_buffer, 3);
//                _field._T._container.assign((Real *)scratch_buffer.data,
//                                            (Real *)scratch_buffer.data + _field._T.size());
//            }
//            _field.calculate_velocities(_grid);
//            // Communicate velocities
//            Communication::communicate(&params, _field.u_matrix());
//            Communication::communicate(&params, _field.v_matrix());
//
//            // Output u,v,p
//            if (t >= output_counter * _output_freq) {
//                output_vtk(timestep, params);
//                output_counter++;
//            }
//            dt = *(Real *)dt_buffer.data;
//            t += dt;
//            timestep++;
//            auto c_end = std::chrono::high_resolution_clock::now();
//            auto time = std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start);
//            // auto time = std::chrono::duration_cast<std::chrono::milliseconds>(ptimecpu);
//            printf(
//                "\rIter: %d, first pass time %.2f ms, pressure time %.2f ms, post pressure %.2f ms, CPU total: %ld
//                ms", it, end - begin, pressure_time, post_end - post_begin, time.count());
//            // Print progress bar
//            logger.progress_bar(t, _t_end);
//        }
//        scratch_buffer.destroy();
//        rs_buffer.destroy();
//        p_buffer.destroy();
//        residual_buffer.destroy();
//        f_buffer.destroy();
//        g_buffer.destroy();
//        u_buffer.destroy();
//        v_buffer.destroy();
//        ubo_buffer.destroy();
//        cell_buffer.destroy();
//        a_data_buffer.destroy();
//        a_offset_buffer.destroy();
//#if ENABLE_PRECOND
//        m_data_buffer.destroy();
//        m_offset_buffer.destroy();
//        z_buffer.destroy();
//#endif
//        d_buffer.destroy();
//        spmv_result_buffer.destroy();
//        counter_buffer.destroy();
//        r_buffer.destroy();
//        deltas_buffer.destroy();
//        neighborhood_buffer.destroy();
//        dt_buffer.destroy();
//        u_boundary_matrix_buffer.destroy();
//        v_boundary_matrix_buffer.destroy();
//        u_rhs_buffer.destroy();
//        v_rhs_buffer.destroy();
//        u_row_start.destroy();
//        v_row_start.destroy();
//        u_col_index.destroy();
//        v_col_index.destroy();
//        if (_calc_temp) {
//            t_new_buffer.destroy();
//            t_old_buffer.destroy();
//            t_rhs_buffer.destroy();
//            t_row_start.destroy();
//            t_col_index.destroy();
//            t_boundary_matrix_buffer.destroy();
//        }
//        // Print Summary
//        if (params.world_rank == 0) logger.finish();
//        // Output u,v,p
//        output_vtk(timestep, params);
//        std::vector<Pipeline> pipelines_to_destroy = {
//            fg_pipeline,
//            rs_pipeline,
//            vel_pipeline,
//            p_pipeline_red,
//            p_pipeline_black,
//            residual_pipeline,
//            p_boundary_pipeline,
//            v_boundary_pipeline,
//            fg_boundary_pipeline,
//            spmv_a_pipeline,
//            saxpy_0_pipeline,
//            saxpy_1_pipeline,
//            saxpy_2_pipeline,
//            reduce_pipeline,
//            vec_dot_vec_0_pipeline,
//            vec_dot_vec_1_pipeline,
//#if ENABLE_PRECOND
//            vec_dot_vec_2_pipeline,
//            spmv_m_pipeline,
//            saxpy_3_pipeline,
//#endif
//            div_pipeline,
//            div_store_pipeline,
//            inc_pipeline,
//            negate_pipeline,
//            min_max_uv_pipeline,
//            reduce_u_pipeline,
//            reduce_v_pipeline,
//            calc_dt_pipeline,
//            boundary_uv_branchless_pipeline
//        };
//        if (_calc_temp) {
//            pipelines_to_destroy.push_back(calc_t_pipeline);
//            pipelines_to_destroy.push_back(boundary_t_branchless);
//        }
//        simulation.cleanup(pipelines_to_destroy);
//    }
//}

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
    VTK_Array *Pressure = VTK_Array::New();
    Pressure->SetName("pressure");
    Pressure->SetNumberOfComponents(1);

    // Velocity Array
    VTK_Array *Velocity = VTK_Array::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);

    // Print pressure and place ghost cells
    for (int j = 1; j < _solver->_grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _solver->_grid.domain().size_x + 1; i++) {
            Real pressure = _solver->_field.p(i, j);
            Pressure->InsertNextTuple(&pressure);

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

    // Add Velocity to Structured Grid
    structuredGrid->GetPointData()->AddArray(Velocity);

    // Add Temperature to Structured Grid
    if (_solver->_calc_temp) {
        VTK_Array *Temperature = VTK_Array::New();
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

void Simulation::build_domain(Domain &domain, int imax_domain, int jmax_domain) {
    domain.imin = 0;
    domain.jmin = 0;
    domain.imax = imax_domain + 2;
    domain.jmax = jmax_domain + 2;
    domain.size_x = imax_domain;
    domain.size_y = jmax_domain;
}
