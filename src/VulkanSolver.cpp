#include "VulkanSolver.hpp"
#define ENABLE_PRECOND 0

void VulkanSolver::solve_pre_pressure(Real &dt) {
    simulation.run_command_buffer(0);
    dt = *(Real *)dt_buffer.data;
}

void VulkanSolver::solve_pressure(Real &res, uint32_t &it) {
    Real delta_new = *(Real *)scratch_buffer.data;
    Real delta_old = delta_new;
    Real delta_zero = delta_new;
    Real cond = _tolerance * _tolerance * delta_zero;
    uint8_t sem_idx = 0;
    VkBufferCopy copy_region;
    copy_region.srcOffset = 0;
    copy_region.dstOffset = 0;
    copy_region.size = deltas_buffer.size;
    it = 0;
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
        simulation.begin_recording(3, 0);
        barrier(simulation, deltas_buffer, 3, VK_ACCESS_SHADER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT);
        vkCmdCopyBuffer(simulation.context.command_buffer[3], deltas_buffer.handle, scratch_buffer.handle, 1,
                        &copy_region);
        simulation.end_recording(3);
        simulation.run_command_buffer(3);
        delta_new = *(Real *)scratch_buffer.data;
    }
}

void VulkanSolver::solve_post_pressure() {
    simulation.run_command_buffer(2);
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
}

void VulkanSolver::initialize() {
    auto grid_x = _grid.imaxb();
    auto grid_y = _grid.jmaxb();
    auto grid_size = grid_x * grid_y;

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

    DiagonalSparseMatrix<Real> A_matrix_diag =
        create_diagonal_matrix(static_cast<PCG *>(_pressure_solver.get())->A, _grid.imaxb(), _grid.jmaxb(),
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
            create_preconditioner_spai(static_cast<PCG *>(_pressure_solver.get())->A, _grid.imaxb(), _grid.jmaxb());
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
    fg_pipeline = simulation.create_compute_pipeline("src/shaders/calc_fg.comp.spv", _calc_temp);
    rs_pipeline = simulation.create_compute_pipeline("src/shaders/calc_rs.comp.spv");
    vel_pipeline = simulation.create_compute_pipeline("src/shaders/calc_vel.comp.spv");
    p_pipeline_red = simulation.create_compute_pipeline("src/shaders/calc_p_redblack_gs.comp.spv", 0);
    p_pipeline_black = simulation.create_compute_pipeline("src/shaders/calc_p_redblack_gs.comp.spv", 1);
    residual_pipeline = simulation.create_compute_pipeline("src/shaders/calc_r.comp.spv");
    p_boundary_pipeline = simulation.create_compute_pipeline("src/shaders/boundary_p.comp.spv");
    v_boundary_pipeline = simulation.create_compute_pipeline("src/shaders/boundary_v.comp.spv");
    fg_boundary_pipeline = simulation.create_compute_pipeline("src/shaders/boundary_fg.comp.spv");
    spmv_a_pipeline = simulation.create_compute_pipeline("src/shaders/spmv.comp.spv", 0);
    saxpy_0_pipeline = simulation.create_compute_pipeline("src/shaders/saxpy.comp.spv", 0);
    saxpy_1_pipeline = simulation.create_compute_pipeline("src/shaders/saxpy.comp.spv", 1);
    saxpy_2_pipeline = simulation.create_compute_pipeline("src/shaders/saxpy.comp.spv", 2);
    vec_dot_vec_0_pipeline = simulation.create_compute_pipeline("src/shaders/vec_dot_vec.comp.spv", 0);
    vec_dot_vec_1_pipeline = simulation.create_compute_pipeline("src/shaders/vec_dot_vec.comp.spv", 1);
#if ENABLE_PRECOND
    vec_dot_vec_2_pipeline = simulation.create_compute_pipeline("src/shaders/vec_dot_vec.comp.spv", 2);
    spmv_m_pipeline = simulation.create_compute_pipeline("src/shaders/spmv.comp.spv", 1);
    saxpy_3_pipeline = simulation.create_compute_pipeline("src/shaders/saxpy.comp.spv", 3);
#endif
    div_pipeline = simulation.create_compute_pipeline("src/shaders/div.comp.spv", 0);
    div_store_pipeline = simulation.create_compute_pipeline("src/shaders/div.comp.spv", 1);
    reduce_pipeline = simulation.create_compute_pipeline("src/shaders/reduce.comp.spv", 0);
    inc_pipeline = simulation.create_compute_pipeline("src/shaders/increment.comp.spv");
    negate_pipeline = simulation.create_compute_pipeline("src/shaders/negate.comp.spv");
    min_max_uv_pipeline = simulation.create_compute_pipeline("src/shaders/min_max_uv.comp.spv");
    reduce_u_pipeline = simulation.create_compute_pipeline("src/shaders/reduce_uv.comp.spv", 0);
    reduce_v_pipeline = simulation.create_compute_pipeline("src/shaders/reduce_uv.comp.spv", 1);
    calc_dt_pipeline = simulation.create_compute_pipeline("src/shaders/calc_dt.comp.spv", _calc_temp);
    boundary_uv_branchless_pipeline = simulation.create_compute_pipeline("src/shaders/boundary_uv_branchless.comp.spv");
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

    /*  dt_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
       VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_SHARING_MODE_EXCLUSIVE, sizeof(Real));*/
    dt_buffer.create(&simulation.context, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
                     VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                     VK_SHARING_MODE_EXCLUSIVE, sizeof(Real));

    auto pcg_solver = (PCG *)_pressure_solver.get();
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
    record_simulation_step(0);
    record_conjugate_gradient_solver(1);
    record_post_pressure(2);
}

VulkanSolver::~VulkanSolver() {
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

void VulkanSolver::record_simulation_step(int command_idx) {
    auto grid_x = _grid.imaxb();
    auto grid_y = _grid.jmaxb();
    auto grid_size = grid_x * grid_y;
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
    // Compute F & G and enforce boundary conditions
    simulation.record_command_buffer(fg_boundary_pipeline, command_idx, 32, 32, grid_x, grid_y);
    if (_calc_temp) {
        t_new_buffer.copy(t_old_buffer, command_idx, false);
        simulation.record_command_buffer(boundary_t_branchless, command_idx, 1024, 1, grid_size, 1);
        barrier(simulation, t_new_buffer, command_idx);
        simulation.record_command_buffer(calc_t_pipeline, command_idx, 32, 32, grid_x, grid_y);
        barrier(simulation, t_new_buffer, command_idx);
    }

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
}

void VulkanSolver::record_conjugate_gradient_solver(int command_idx) {
    auto grid_x = _grid.imaxb();
    auto grid_y = _grid.jmaxb();
    auto grid_size = grid_x * grid_y;
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
    // barrier(simulation, p_buffer, command_idx);
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
}

void VulkanSolver::record_post_pressure(int command_idx) {
    auto grid_x = _grid.imaxb();
    auto grid_y = _grid.jmaxb();
    auto grid_size = grid_x * grid_y;
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
}
