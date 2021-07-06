#include "VulkanUtils.hpp"

void scalar_div(GPUSimulation &simulation, Pipeline &pipeline, int command_idx) {
    simulation.record_command_buffer(pipeline, command_idx, 1, 1, 1, 1);
}

void vec_saxpy(GPUSimulation &simulation, Pipeline &pipeline, int command_idx, int dim) {
    simulation.record_command_buffer(pipeline, command_idx, 1024, 1, dim, 1);
}

void barrier(GPUSimulation &simulation, Buffer &buffer, int command_idx, VkAccessFlags src_access,
             VkAccessFlags dst_access) {
    VkBufferMemoryBarrier res_barrier = buffer_barrier(buffer.handle, src_access, dst_access);
    vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                         VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &res_barrier, 0, 0);
}

Real vec_dp_immediate(GPUSimulation &simulation, Buffer &v1, Buffer &v2, Buffer &residual_buffer,
                      Buffer &scratch_buffer, Pipeline &pipeline, Pipeline &reduce_pipeline, int command_idx, int dim) {
    simulation.begin_recording(command_idx);
    vkCmdFillBuffer(simulation.context.command_buffer[command_idx], residual_buffer.handle, 0, residual_buffer.size, 0);
    VkBufferMemoryBarrier fill_barrier =
        buffer_barrier(residual_buffer.handle, VK_ACCESS_TRANSFER_WRITE_BIT, VK_ACCESS_SHADER_WRITE_BIT);
    vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_TRANSFER_BIT,
                         VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &fill_barrier, 0, 0);
    simulation.record_command_buffer(pipeline, command_idx, 1024, 1, dim, 1);
    VkBufferMemoryBarrier res_barrier = buffer_barrier(residual_buffer.handle, VK_ACCESS_SHADER_WRITE_BIT,
                                                       VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_SHADER_READ_BIT);
    vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                         VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &res_barrier, 0, 0);
    simulation.record_command_buffer(reduce_pipeline, command_idx, 32, 1, 32, 1); // TODO
    simulation.end_recording(command_idx);
    simulation.run_command_buffer(command_idx);
    residual_buffer.copy(scratch_buffer, command_idx);
    return *(Real *)scratch_buffer.data;
}

void vec_dp(GPUSimulation &simulation, Buffer &residual_buffer, Buffer &counter_buffer, Pipeline &pipeline,
            Pipeline &reduce_pipeline, int command_idx, int dim) {
    /* vkCmdFillBuffer(simulation.context.command_buffer[command_idx], residual_buffer.handle, 0, residual_buffer.size,
     0); VkBufferMemoryBarrier fill_barrier = buffer_barrier(residual_buffer.handle, VK_ACCESS_TRANSFER_WRITE_BIT,
     VK_ACCESS_SHADER_WRITE_BIT); vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx],
     VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &fill_barrier, 0, 0);*/
    vkCmdFillBuffer(simulation.context.command_buffer[command_idx], counter_buffer.handle, 0, counter_buffer.size, 1);
    simulation.record_command_buffer(pipeline, command_idx, 1024, 1, dim, 1);
    VkBufferMemoryBarrier res_barrier = buffer_barrier(residual_buffer.handle, VK_ACCESS_SHADER_WRITE_BIT,
                                                       VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_SHADER_READ_BIT);
    vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                         VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &res_barrier, 0, 0);
    int num_wgs = ceil(dim / 1024.0f);

    VkBufferMemoryBarrier fill_barrier = buffer_barrier(counter_buffer.handle, VK_ACCESS_TRANSFER_WRITE_BIT,
                                                        VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_SHADER_READ_BIT);
    VkBufferMemoryBarrier counter_barrier = buffer_barrier(counter_buffer.handle, VK_ACCESS_SHADER_WRITE_BIT,
                                                           VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_SHADER_READ_BIT);
    vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_TRANSFER_BIT,
                         VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &fill_barrier, 0, 0);
    std::vector<VkBufferMemoryBarrier> barriers{res_barrier, fill_barrier};
    while (num_wgs != 1) {
        simulation.record_command_buffer(reduce_pipeline, command_idx, 1024, 1, num_wgs, 1);
        num_wgs = ceil(num_wgs / 1024.0f);
        if (num_wgs > 1) {
            vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                                 VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 2, barriers.data(), 0, 0);
        }
    }
}

void calc_residual(GPUSimulation &simulation, Buffer &residual_buffer, Buffer &counter_buffer,
                   Pipeline &residual_pipeline, Pipeline &reduce_pipeline, int command_idx, int grid_x, int grid_y) {
    vkCmdFillBuffer(simulation.context.command_buffer[command_idx], residual_buffer.handle, 0, residual_buffer.size, 0);
    vkCmdFillBuffer(simulation.context.command_buffer[command_idx], counter_buffer.handle, 0, counter_buffer.size, 0);
    VkBufferMemoryBarrier fill_res_barrier = buffer_barrier(residual_buffer.handle, VK_ACCESS_TRANSFER_WRITE_BIT,
                                                            VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_SHADER_READ_BIT);
    vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_TRANSFER_BIT,
                         VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &fill_res_barrier, 0, 0);
    simulation.record_command_buffer(residual_pipeline, command_idx, 32, 32, grid_x, grid_y);
    VkBufferMemoryBarrier res_barrier = buffer_barrier(residual_buffer.handle, VK_ACCESS_SHADER_WRITE_BIT,
                                                       VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_SHADER_READ_BIT);
    vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                         VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &res_barrier, 0, 0);
    auto grid_dim = grid_x * grid_y;
    int num_wgs = grid_dim;
    VkBufferMemoryBarrier fill_barrier = buffer_barrier(counter_buffer.handle, VK_ACCESS_TRANSFER_WRITE_BIT,
                                                        VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_SHADER_READ_BIT);

    VkBufferMemoryBarrier counter_barrier = buffer_barrier(counter_buffer.handle, VK_ACCESS_SHADER_WRITE_BIT,
                                                           VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_SHADER_READ_BIT);
    vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_TRANSFER_BIT,
                         VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &fill_barrier, 0, 0);
    std::vector<VkBufferMemoryBarrier> barriers{res_barrier, counter_barrier};
    while (num_wgs != 1) {
        simulation.record_command_buffer(reduce_pipeline, command_idx, 1024, 1, num_wgs, 1);
        num_wgs = ceil(num_wgs / 1024.0f);
        if (num_wgs > 1) {
            vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                                 VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 2, barriers.data(), 0, 0);
        }
    }
}

void uv_max(GPUSimulation &simulation, Buffer &residual_buffer, Buffer &counter_buffer, Pipeline &min_max_uv_pipeline,
            Pipeline &reduce_u_pipeline, Pipeline &reduce_v_pipeline, int command_idx, int dim) {
    vkCmdFillBuffer(simulation.context.command_buffer[command_idx], residual_buffer.handle, 0, residual_buffer.size, 0);
    VkBufferMemoryBarrier fill_barrier =
        buffer_barrier(residual_buffer.handle, VK_ACCESS_TRANSFER_WRITE_BIT, VK_ACCESS_SHADER_WRITE_BIT);
    vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_TRANSFER_BIT,
                         VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &fill_barrier, 0, 0);
    VkBufferMemoryBarrier res_barrier = buffer_barrier(residual_buffer.handle, VK_ACCESS_SHADER_WRITE_BIT,
                                                       VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_SHADER_READ_BIT);
    simulation.record_command_buffer(min_max_uv_pipeline, command_idx, 1024, 1, dim, 1);
    vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                         VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &res_barrier, 0, 0);
    int num_wgs = ceil(dim / 1024.0f);
    vkCmdFillBuffer(simulation.context.command_buffer[command_idx], counter_buffer.handle, 0, counter_buffer.size, 1);
    fill_barrier = buffer_barrier(counter_buffer.handle, VK_ACCESS_TRANSFER_WRITE_BIT,
                                  VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_SHADER_READ_BIT);
    VkBufferMemoryBarrier counter_barrier = buffer_barrier(counter_buffer.handle, VK_ACCESS_SHADER_WRITE_BIT,
                                                           VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_SHADER_READ_BIT);
    vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_TRANSFER_BIT,
                         VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &fill_barrier, 0, 0);
    std::vector<VkBufferMemoryBarrier> barriers{res_barrier, counter_barrier};
    while (num_wgs != 1) {
        simulation.record_command_buffer(reduce_u_pipeline, command_idx, 1024, 1, num_wgs, 1);
        simulation.record_command_buffer(reduce_v_pipeline, command_idx, 1024, 1, num_wgs, 1);
        num_wgs = ceil(num_wgs / 1024.0f);
        if (num_wgs > 1) {
            vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                                 VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 2, barriers.data(), 0, 0);
        }
    }
}

VkDescriptorSetLayoutBinding descriptor_set_layout_binding(VkDescriptorType type, VkShaderStageFlags stageFlags,
                                                           uint32_t binding, uint32_t descriptorCount) {
    VkDescriptorSetLayoutBinding setLayoutBinding{};
    setLayoutBinding.descriptorType = type;
    setLayoutBinding.stageFlags = stageFlags;
    setLayoutBinding.binding = binding;
    setLayoutBinding.descriptorCount = descriptorCount;
    return setLayoutBinding;
}

VkWriteDescriptorSet write_descriptor_set(VkDescriptorSet dstSet, VkDescriptorType type, uint32_t binding,
                                          VkDescriptorBufferInfo *bufferInfo, uint32_t descriptorCount) {
    VkWriteDescriptorSet write_descriptor_set{};
    write_descriptor_set.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
    write_descriptor_set.dstSet = dstSet;
    write_descriptor_set.descriptorType = type;
    write_descriptor_set.dstBinding = binding;
    write_descriptor_set.pBufferInfo = bufferInfo;
    write_descriptor_set.descriptorCount = descriptorCount;
    return write_descriptor_set;
}

uint32_t find_memory_type(VkPhysicalDevice physical_device, VkMemoryPropertyFlags properties,
                          uint32_t memory_type_bits) {
    VkPhysicalDeviceMemoryProperties memory_properties;
    vkGetPhysicalDeviceMemoryProperties(physical_device, &memory_properties);
    for (uint32_t i = 0; i < memory_properties.memoryTypeCount; ++i) {
        if ((memory_type_bits & (1 << i)) &&
            ((memory_properties.memoryTypes[i].propertyFlags & properties) == properties))
            return i;
    }
    return -1;
}

VkBufferMemoryBarrier buffer_barrier(VkBuffer handle, VkAccessFlags src_access_mask, VkAccessFlags dst_access_mask) {
    VkBufferMemoryBarrier result = {VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER};
    result.srcAccessMask = src_access_mask;
    result.dstAccessMask = dst_access_mask;
    result.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    result.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    result.buffer = handle;
    result.offset = 0;
    result.size = VK_WHOLE_SIZE;
    return result;
}