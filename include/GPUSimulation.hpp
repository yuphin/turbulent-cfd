#include <vector>
#include <vulkan/vulkan.h>
#ifdef NDEBUG
constexpr bool enable_validation = false;
#else
constexpr bool enable_validation = true;
#endif

enum CommandBufferIndex {
    COMMAND_BUFFER_IDX_0,
    COMMAND_BUFFER_IDX_1,
    COMMAND_BUFFER_IDX_2,
    COMMAND_BUFFER_IDX_3,
};

constexpr int MAX_DESCRIPTOR_SETS = 4;

struct UBOData {
    int imax;
    int jmax;
    int size;
    float nu;
    float alpha;
    float beta;
    float gx;
    float gy;
    float dx;
    float dy;
    float dx2;
    float dy2;
    float inv_dx;
    float inv_dy;
    float gamma;
    float PI;
    float UI;
    float VI;
    float tau;
    int num_diags;
};

#define VK_CHECK(f)                                                                                                    \
    {                                                                                                                  \
        VkResult res = (f);                                                                                            \
        if (res != VK_SUCCESS) {                                                                                       \
            printf("Fatal : VkResult is %d in %s at line %d\n", res, __FILE__, __LINE__);                              \
            assert(res == VK_SUCCESS);                                                                                 \
        }                                                                                                              \
    }

static VKAPI_ATTR VkBool32 VKAPI_CALL debugReportCallbackFn(VkDebugReportFlagsEXT flags,
                                                            VkDebugReportObjectTypeEXT objectType, uint64_t object,
                                                            size_t location, int32_t messageCode,
                                                            const char *pLayerPrefix, const char *pMessage,
                                                            void *pUserData) {
    // Work around errors regarding specialization constants and different buffers in the same shader
    if (strstr(pMessage, "The Vulkan spec states: Descriptors in each bound descriptor set, specified via "
                         "vkCmdBindDescriptorSets, must be valid"))
        return VK_FALSE;
    // Work around errors regarding unupdated descriptor sets
    if (strstr(pMessage, "is being used in draw but has not been updated")) 
        return VK_FALSE;
    printf("Debug Report: %s: %s\n", pLayerPrefix, pMessage);
    return VK_FALSE;
}

VkDescriptorSetLayoutBinding descriptor_set_layout_binding(VkDescriptorType type, VkShaderStageFlags stageFlags,
                                                           uint32_t binding, uint32_t descriptorCount = 1) {
    VkDescriptorSetLayoutBinding setLayoutBinding{};
    setLayoutBinding.descriptorType = type;
    setLayoutBinding.stageFlags = stageFlags;
    setLayoutBinding.binding = binding;
    setLayoutBinding.descriptorCount = descriptorCount;
    return setLayoutBinding;
}

VkWriteDescriptorSet write_descriptor_set(VkDescriptorSet dstSet, VkDescriptorType type, uint32_t binding,
                                          VkDescriptorBufferInfo *bufferInfo, uint32_t descriptorCount = 1) {
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

struct Pipeline {
    VkPipeline pipeline;
    VkPipelineLayout pipeline_layout;
    VkShaderModule shader_module;
};
struct Context {
    VkDevice device;
    VkQueue compute_queue;
    VkPhysicalDevice physical_device;
    VkCommandPool command_pool;
    std::vector<VkCommandBuffer> command_buffer;
};

struct Buffer {
    Context *ctx = nullptr;
    VkBuffer handle = 0;
    VkDeviceMemory memory = 0;
    void *data = nullptr;
    size_t size;
    VkDeviceSize alignment = 0;
    VkBufferUsageFlags usage_flags = 0;
    VkMemoryPropertyFlags mem_prop_flags = 0;
    VkDescriptorBufferInfo descriptor = {};
    void destroy() {
        if (handle) vkDestroyBuffer(ctx->device, handle, nullptr);
        if (memory) vkFreeMemory(ctx->device, memory, nullptr);
    }
    void create(Context *ctx, VkBufferUsageFlags usage, VkMemoryPropertyFlags mem_property_flags,
                VkSharingMode sharing_mode, VkDeviceSize size, void *data = nullptr, bool use_staging = false,
                int command_idx = 0) {
        if (!this->ctx) {
            this->ctx = ctx;
            this->mem_prop_flags = mem_property_flags;
            this->usage_flags = usage;
        }
        if (use_staging) {
            Buffer staging_buffer;
            assert(mem_property_flags & VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
            staging_buffer.create(ctx, VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
                                  VK_MEMORY_PROPERTY_HOST_COHERENT_BIT | VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT,
                                  VK_SHARING_MODE_EXCLUSIVE, size, data, command_idx);
            this->create(ctx, VK_BUFFER_USAGE_TRANSFER_DST_BIT | usage, mem_property_flags, sharing_mode, size, nullptr,
                         /* use_staging */ false, command_idx);

            // @performance
            // VK_CHECK(vkResetCommandPool(ctx->device, ctx->command_pool, 0));
            VkCommandBufferBeginInfo begin_info = {};
            begin_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
            begin_info.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
            VK_CHECK(vkBeginCommandBuffer(ctx->command_buffer[command_idx], &begin_info));
            VkBufferCopy copy_region = {0, 0, VkDeviceSize(size)};
            vkCmdCopyBuffer(ctx->command_buffer[command_idx], staging_buffer.handle, this->handle, 1, &copy_region);
            VkBufferMemoryBarrier copy_barrier =
                buffer_barrier(handle, VK_ACCESS_TRANSFER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT);
            vkCmdPipelineBarrier(ctx->command_buffer[command_idx], VK_PIPELINE_STAGE_TRANSFER_BIT,
                                 VK_PIPELINE_STAGE_ALL_COMMANDS_BIT, VK_DEPENDENCY_BY_REGION_BIT, 0, 0, 1,
                                 &copy_barrier, 0, 0);
            vkEndCommandBuffer(ctx->command_buffer[command_idx]);
            VkSubmitInfo submitInfo = {VK_STRUCTURE_TYPE_SUBMIT_INFO};
            submitInfo.commandBufferCount = 1;
            submitInfo.pCommandBuffers = &ctx->command_buffer[command_idx];

            VK_CHECK(vkQueueSubmit(ctx->compute_queue, 1, &submitInfo, VK_NULL_HANDLE));

            VK_CHECK(vkDeviceWaitIdle(ctx->device));
            staging_buffer.destroy();
        } else {
            // Create the buffer handle
            VkBufferCreateInfo buffer_create_info = {VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
            buffer_create_info.usage = usage;
            buffer_create_info.size = size;
            buffer_create_info.sharingMode = sharing_mode;
            VK_CHECK(vkCreateBuffer(ctx->device, &buffer_create_info, nullptr, &handle));

            // Create the memory backing up the buffer handle
            VkMemoryRequirements mem_reqs;
            VkMemoryAllocateInfo mem_alloc_info = {VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
            vkGetBufferMemoryRequirements(ctx->device, this->handle, &mem_reqs);

            mem_alloc_info.allocationSize = mem_reqs.size;
            // Find a memory type index that fits the properties of the buffer
            mem_alloc_info.memoryTypeIndex =
                find_memory_type(ctx->physical_device, mem_property_flags, mem_reqs.memoryTypeBits);
            VK_CHECK(vkAllocateMemory(ctx->device, &mem_alloc_info, nullptr, &memory));

            alignment = mem_reqs.alignment;
            this->size = size;
            usage_flags = usage;
            if (mem_property_flags & VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT) {
                map();
            }
            // If a pointer to the buffer data has been passed, map the buffer and copy over the data
            if (data != nullptr) {
                memcpy(this->data, data, size);
                if ((mem_property_flags & VK_MEMORY_PROPERTY_HOST_COHERENT_BIT) == 0) {
                    flush();
                }
            }
            // Initialize a default descriptor that covers the whole buffer size
            this->prepare_descriptor();
            this->bind();
        }
    }
    void prepare_descriptor(VkDeviceSize size = VK_WHOLE_SIZE, VkDeviceSize offset = 0) {
        descriptor.offset = offset;
        descriptor.buffer = handle;
        descriptor.range = size;
    }
    void bind(VkDeviceSize offset = 0) { VK_CHECK(vkBindBufferMemory(ctx->device, handle, memory, offset)); }
    void map() { vkMapMemory(ctx->device, memory, 0, VK_WHOLE_SIZE, 0, &data); }
    void unmap() { vkUnmapMemory(ctx->device, memory); }
    void flush(VkDeviceSize size = VK_WHOLE_SIZE, VkDeviceSize offset = 0) {
        VkMappedMemoryRange mapped_range = {};
        mapped_range.sType = VK_STRUCTURE_TYPE_MAPPED_MEMORY_RANGE;
        mapped_range.memory = memory;
        mapped_range.offset = offset;
        mapped_range.size = size;
        VK_CHECK(vkFlushMappedMemoryRanges(ctx->device, 1, &mapped_range));
    }
    void upload(VkDeviceSize size, int command_idx, void *data = nullptr, Buffer *staging_buffer = nullptr) {
        if ((mem_prop_flags & (VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT))) {
            memcpy(this->data, data, size);
        } else {
            // Assume device local memory
            staging_buffer->upload(size, command_idx, data);
            // VK_CHECK(vkResetCommandPool(ctx->device, ctx->command_pool, 0));
            VkCommandBufferBeginInfo begin_info = {};
            begin_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
            begin_info.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
            VK_CHECK(vkBeginCommandBuffer(ctx->command_buffer[command_idx], &begin_info));
            VkBufferCopy copy_region = {0, 0, VkDeviceSize(size)};
            vkCmdCopyBuffer(ctx->command_buffer[command_idx], staging_buffer->handle, this->handle, 1, &copy_region);
            VkBufferMemoryBarrier copy_barrier =
                buffer_barrier(handle, VK_ACCESS_TRANSFER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT);
            vkCmdPipelineBarrier(ctx->command_buffer[command_idx], VK_PIPELINE_STAGE_TRANSFER_BIT,
                                 VK_PIPELINE_STAGE_ALL_COMMANDS_BIT, VK_DEPENDENCY_BY_REGION_BIT, 0, 0, 1,
                                 &copy_barrier, 0, 0);
            vkEndCommandBuffer(ctx->command_buffer[command_idx]);
            VkSubmitInfo submitInfo = {VK_STRUCTURE_TYPE_SUBMIT_INFO};
            submitInfo.commandBufferCount = 1;
            submitInfo.pCommandBuffers = &ctx->command_buffer[command_idx];
            VK_CHECK(vkQueueSubmit(ctx->compute_queue, 1, &submitInfo, VK_NULL_HANDLE));
            VK_CHECK(vkDeviceWaitIdle(ctx->device));
        }
    }
    void copy(Buffer &dst_buffer, int command_idx, bool reset = true, VkDeviceSize src_offset = 0,
              VkDeviceSize dst_offset = 0, VkDeviceSize copy_size = 0) {
        VkBufferCopy copy_region;
        copy_region.srcOffset = 0;
        copy_region.dstOffset = 0;
        copy_region.size = copy_size == 0 ? size : copy_size;
        VkCommandBufferBeginInfo begin_info = {};
        begin_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
        begin_info.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
        if (reset) {
            // VK_CHECK(vkResetCommandPool(ctx->device, ctx->command_pool, 0));
            VK_CHECK(vkBeginCommandBuffer(ctx->command_buffer[command_idx], &begin_info));
        }

        vkCmdCopyBuffer(ctx->command_buffer[command_idx], handle, dst_buffer.handle, 1, &copy_region);
        VkBufferMemoryBarrier copy_barrier =
            buffer_barrier(handle, VK_ACCESS_TRANSFER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT);
        vkCmdPipelineBarrier(ctx->command_buffer[command_idx], VK_PIPELINE_STAGE_TRANSFER_BIT,
                             VK_PIPELINE_STAGE_ALL_COMMANDS_BIT, VK_DEPENDENCY_BY_REGION_BIT, 0, 0, 1, &copy_barrier, 0,
                             0);
        if (reset) {
            vkEndCommandBuffer(ctx->command_buffer[command_idx]);
            VkSubmitInfo submitInfo = {VK_STRUCTURE_TYPE_SUBMIT_INFO};
            submitInfo.commandBufferCount = 1;
            submitInfo.pCommandBuffers = &ctx->command_buffer[command_idx];
            VK_CHECK(vkQueueSubmit(ctx->compute_queue, 1, &submitInfo, VK_NULL_HANDLE));
            VK_CHECK(vkDeviceWaitIdle(ctx->device));
        }
    }
};

struct Descriptor {
    Buffer handle;
    int binding;
    int set_idx = 0;
    VkDescriptorType desc_type = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
};

class GPUSimulation {
  public:
    void init(UBOData &data) {
        create_instance();
        find_physical_device();
        create_device();
        create_descriptor_pool();
        this->data = data;
        descriptor_sets.resize(MAX_DESCRIPTOR_SETS);
        descriptor_set_layouts.resize(MAX_DESCRIPTOR_SETS);
    }
    void create_instance() {
        std::vector<const char *> enabled_extensions;

        if (enable_validation) {
            uint32_t layer_count;
            vkEnumerateInstanceLayerProperties(&layer_count, NULL);

            std::vector<VkLayerProperties> layer_props(layer_count);
            vkEnumerateInstanceLayerProperties(&layer_count, layer_props.data());

            bool found_layer = false;
            for (VkLayerProperties prop : layer_props) {

                if (strcmp("VK_LAYER_KHRONOS_validation", prop.layerName) == 0) {
                    found_layer = true;
                    break;
                }
            }

            if (!found_layer) {
                throw std::runtime_error("Layer  VK_LAYER_KHRONOS_validation not supported\n");
            }
            enabled_layers.push_back("VK_LAYER_KHRONOS_validation"); // Alright, we can use this layer.

            uint32_t extension_count;
            vkEnumerateInstanceExtensionProperties(NULL, &extension_count, NULL);
            std::vector<VkExtensionProperties> extension_props(extension_count);
            vkEnumerateInstanceExtensionProperties(NULL, &extension_count, extension_props.data());

            bool found_extension = false;
            for (VkExtensionProperties prop : extension_props) {
                if (strcmp(VK_EXT_DEBUG_REPORT_EXTENSION_NAME, prop.extensionName) == 0) {
                    found_extension = true;
                    break;
                }
            }

            if (!found_extension) {
                throw std::runtime_error("Extension VK_EXT_DEBUG_REPORT_EXTENSION_NAME not supported\n");
            }
            enabled_extensions.push_back(VK_EXT_DEBUG_REPORT_EXTENSION_NAME);
        }

        VkApplicationInfo application_info = {};
        application_info.sType = VK_STRUCTURE_TYPE_APPLICATION_INFO;
        application_info.pApplicationName = "Fluidchen";
        application_info.applicationVersion = 0;
        application_info.pEngineName = "Fluidchen";
        application_info.engineVersion = 0;
        application_info.apiVersion = VK_API_VERSION_1_2;

        VkInstanceCreateInfo create_info = {};
        create_info.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
        create_info.flags = 0;
        create_info.pApplicationInfo = &application_info;

        create_info.enabledLayerCount = enabled_layers.size();
        create_info.ppEnabledLayerNames = enabled_layers.data();
        create_info.enabledExtensionCount = enabled_extensions.size();
        create_info.ppEnabledExtensionNames = enabled_extensions.data();

        VK_CHECK(vkCreateInstance(&create_info, NULL, &instance));

        if (enable_validation) {
            VkDebugReportCallbackCreateInfoEXT create_info = {};
            create_info.sType = VK_STRUCTURE_TYPE_DEBUG_REPORT_CALLBACK_CREATE_INFO_EXT;
            create_info.flags = VK_DEBUG_REPORT_ERROR_BIT_EXT | VK_DEBUG_REPORT_WARNING_BIT_EXT |
                                VK_DEBUG_REPORT_PERFORMANCE_WARNING_BIT_EXT;
            create_info.pfnCallback = &debugReportCallbackFn;

            auto vkCreateDebugReportCallbackEXT =
                (PFN_vkCreateDebugReportCallbackEXT)vkGetInstanceProcAddr(instance, "vkCreateDebugReportCallbackEXT");
            if (vkCreateDebugReportCallbackEXT == nullptr) {
                throw std::runtime_error("Could not load vkCreateDebugReportCallbackEXT");
            }

            VK_CHECK(vkCreateDebugReportCallbackEXT(instance, &create_info, NULL, &debug_report_callback));
        }
    }

    void find_physical_device() {
        uint32_t device_count;
        vkEnumeratePhysicalDevices(instance, &device_count, NULL);
        if (device_count == 0) {
            throw std::runtime_error("could not find a device with vulkan support");
        }

        std::vector<VkPhysicalDevice> devices(device_count);
        vkEnumeratePhysicalDevices(instance, &device_count, devices.data());

        // TODO: Check device limitations
        for (VkPhysicalDevice device : devices) {
            if (true) { // As above stated, we do no feature checks, so just accept.
                context.physical_device = device;
                break;
            }
        }
    }

    uint32_t get_compute_queue_family_index() {

        uint32_t queue_family_count;
        vkGetPhysicalDeviceQueueFamilyProperties(context.physical_device, &queue_family_count, NULL);

        std::vector<VkQueueFamilyProperties> queue_families(queue_family_count);
        vkGetPhysicalDeviceQueueFamilyProperties(context.physical_device, &queue_family_count, queue_families.data());

        uint32_t i = 0;
        for (; i < queue_families.size(); ++i) {
            VkQueueFamilyProperties props = queue_families[i];

            if (props.queueCount > 0 && (props.queueFlags & VK_QUEUE_COMPUTE_BIT)) {
                // Found
                break;
            }
        }

        if (i == queue_families.size()) {
            throw std::runtime_error("Could not find a queue family that supports compute operations");
        }

        return i;
    }

    void create_device() {
        VkDeviceQueueCreateInfo queue_create_info = {};
        queue_create_info.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
        compute_queue_family_index = get_compute_queue_family_index();
        queue_create_info.queueFamilyIndex = compute_queue_family_index;
        queue_create_info.queueCount = 1;
        float queue_priorities = 1.0;
        queue_create_info.pQueuePriorities = &queue_priorities;

        VkDeviceCreateInfo device_create_info = {};

        // None for now
        VkPhysicalDeviceFeatures device_features = {};

        device_create_info.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
        device_create_info.enabledLayerCount = enabled_layers.size();
        device_create_info.ppEnabledLayerNames = enabled_layers.data();
        device_create_info.pQueueCreateInfos = &queue_create_info;
        device_create_info.queueCreateInfoCount = 1;
        device_create_info.pEnabledFeatures = &device_features;

        VK_CHECK(vkCreateDevice(context.physical_device, &device_create_info, NULL, &context.device));

        vkGetDeviceQueue(context.device, compute_queue_family_index, 0, &context.compute_queue);
    }
    void create_mem(VkDeviceMemory &mem, VkMemoryPropertyFlags mem_flags, size_t size) {
        uint32_t mem_type_index = find_memory_type(context.physical_device, mem_flags, size);
        assert(mem_type_index != ~0u);
        VkMemoryAllocateInfo allocate_info = {VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
        allocate_info.allocationSize = size;
        allocate_info.memoryTypeIndex = mem_type_index;
        VK_CHECK(vkAllocateMemory(context.device, &allocate_info, 0, &mem));
    }

    void create_descriptor_pool() {
        std::vector<VkDescriptorPoolSize> pool_sizes = {{VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 32},
                                                        {VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 32}};
        VkDescriptorPoolCreateInfo descriptor_pool_create_info = {};
        descriptor_pool_create_info.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
        descriptor_pool_create_info.maxSets = 32;
        descriptor_pool_create_info.poolSizeCount = pool_sizes.size();
        descriptor_pool_create_info.pPoolSizes = pool_sizes.data();

        VK_CHECK(vkCreateDescriptorPool(context.device, &descriptor_pool_create_info, NULL, &descriptor_pool));
    }

    void create_descriptor_set_layout(const std::vector<VkDescriptorSetLayoutBinding> &bindings) {
        VkDescriptorSetLayoutCreateInfo dsl_create_info = {VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
        dsl_create_info.bindingCount = bindings.size();
        dsl_create_info.pBindings = bindings.data();
        for (int i = 0; i < MAX_DESCRIPTOR_SETS; i++) {
            VK_CHECK(vkCreateDescriptorSetLayout(context.device, &dsl_create_info, NULL, &descriptor_set_layouts[i]));
            VkDescriptorSetAllocateInfo descriptor_set_allocate_info = {VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
            descriptor_set_allocate_info.descriptorPool = descriptor_pool;
            descriptor_set_allocate_info.descriptorSetCount = 1;
            descriptor_set_allocate_info.pSetLayouts = &descriptor_set_layouts[i];

            VK_CHECK(vkAllocateDescriptorSets(context.device, &descriptor_set_allocate_info, &descriptor_sets[i]));
        }
    }

    void push_descriptors(const std::vector<Descriptor> &descriptors) {
        std::vector<VkWriteDescriptorSet> write_descriptor_sets;
        for (const auto &descriptor : descriptors) {
            VkWriteDescriptorSet write_descriptor_set = {};
            write_descriptor_set.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            write_descriptor_set.dstSet = descriptor_sets[descriptor.set_idx];
            write_descriptor_set.dstBinding = descriptor.binding;
            write_descriptor_set.descriptorCount = 1;
            write_descriptor_set.descriptorType = descriptor.desc_type;
            write_descriptor_set.pBufferInfo = &descriptor.handle.descriptor;
            write_descriptor_sets.push_back(write_descriptor_set);
        }
        vkUpdateDescriptorSets(context.device, write_descriptor_sets.size(), write_descriptor_sets.data(), 0, NULL);

        /*  vkCmdBindDescriptorSets(command_buffer, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline_layout, 0, 1,
           &descriptor_set, 0, NULL);*/
    }

    uint32_t *read_file(uint32_t &length, const char *filename) {

        FILE *fp = fopen(filename, "rb");
        if (fp == NULL) {
            printf("Could not find or open file: %s\n", filename);
        }

        // get file size.
        fseek(fp, 0, SEEK_END);
        long filesize = ftell(fp);
        fseek(fp, 0, SEEK_SET);

        long filesizepadded = long(ceil(filesize / 4.0)) * 4;

        // read file contents.
        char *str = new char[filesizepadded];
        fread(str, filesize, sizeof(char), fp);
        fclose(fp);

        // data padding.
        for (int i = filesize; i < filesizepadded; i++) {
            str[i] = 0;
        }

        length = filesizepadded;
        return (uint32_t *)str;
    }

    Pipeline create_compute_pipeline(const char *shader_path, uint32_t specialization_data = -1) {
        uint32_t filelength;
        uint32_t *code = read_file(filelength, shader_path);
        VkShaderModuleCreateInfo create_info = {};
        create_info.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
        create_info.pCode = code;
        create_info.codeSize = filelength;
        VkShaderModule compute_shader_module;
        VkPipelineLayout pipeline_layout;
        VkPipeline pipeline;
        VK_CHECK(vkCreateShaderModule(context.device, &create_info, NULL, &compute_shader_module));
        delete[] code;

        VkSpecializationMapEntry entry;
        entry.constantID = 0;
        entry.size = sizeof(uint32_t);
        entry.offset = 0;
        VkSpecializationInfo specialization_info{};
        specialization_info.dataSize = sizeof(uint32_t);
        specialization_info.mapEntryCount = 1;
        specialization_info.pMapEntries = &entry;
        specialization_info.pData = &specialization_data;

        VkPushConstantRange push_constant_range{};
        push_constant_range.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        push_constant_range.offset = 0;
        push_constant_range.size = sizeof(UBOData);

        VkPipelineShaderStageCreateInfo shader_stage_create_info = {};
        shader_stage_create_info.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        shader_stage_create_info.stage = VK_SHADER_STAGE_COMPUTE_BIT;
        shader_stage_create_info.module = compute_shader_module;
        if (specialization_data != -1) {
            shader_stage_create_info.pSpecializationInfo = &specialization_info;
        }
        shader_stage_create_info.pName = "main";

        VkPipelineLayoutCreateInfo pipeline_layout_create_info = {};
        pipeline_layout_create_info.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
        pipeline_layout_create_info.setLayoutCount = MAX_DESCRIPTOR_SETS;
        pipeline_layout_create_info.pSetLayouts = descriptor_set_layouts.data();
        pipeline_layout_create_info.pushConstantRangeCount = 1;
        pipeline_layout_create_info.pPushConstantRanges = &push_constant_range;
        VK_CHECK(vkCreatePipelineLayout(context.device, &pipeline_layout_create_info, NULL, &pipeline_layout));

        VkComputePipelineCreateInfo pipeline_create_info = {};
        pipeline_create_info.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
        pipeline_create_info.stage = shader_stage_create_info;
        pipeline_create_info.layout = pipeline_layout;
        VK_CHECK(vkCreateComputePipelines(context.device, VK_NULL_HANDLE, 1, &pipeline_create_info, NULL, &pipeline));
        return {pipeline, pipeline_layout, compute_shader_module};
    }

    void create_command_pool() {

        VkCommandPoolCreateInfo command_pool_create_info = {};
        command_pool_create_info.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
        command_pool_create_info.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
        command_pool_create_info.queueFamilyIndex = compute_queue_family_index;
        VK_CHECK(vkCreateCommandPool(context.device, &command_pool_create_info, NULL, &context.command_pool));
        VkCommandBufferAllocateInfo command_buffer_allocate_info = {};
        command_buffer_allocate_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
        command_buffer_allocate_info.commandPool = context.command_pool;
        command_buffer_allocate_info.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
        context.command_buffer.resize(4);
        command_buffer_allocate_info.commandBufferCount = 4;
        VK_CHECK(
            vkAllocateCommandBuffers(context.device, &command_buffer_allocate_info, context.command_buffer.data()));
    }

    void begin_recording(int command_idx,
                         VkCommandBufferUsageFlags flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT) {

        // VK_CHECK(vkResetCommandPool(context.device, context.command_pool, 0));
        VkCommandBufferBeginInfo begin_info = {};
        begin_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
        begin_info.flags = flags;
        VK_CHECK(vkBeginCommandBuffer(context.command_buffer[command_idx], &begin_info)); // start recording commands.
    }

    void end_recording(int command_idx) { VK_CHECK(vkEndCommandBuffer(context.command_buffer[command_idx])); }

    void begin_end_record_command_buffer(Pipeline pipeline, int command_idx, int wg_x, int wg_y, int width,
                                         int height) {
        begin_recording(command_idx);
        record_command_buffer(pipeline, command_idx, wg_x, wg_y, width, height);
        end_recording(command_idx);
    }

    void record_command_buffer(Pipeline pipeline, int command_idx, int wg_x, int wg_y, int width, int height) {

        vkCmdBindPipeline(context.command_buffer[command_idx], VK_PIPELINE_BIND_POINT_COMPUTE, pipeline.pipeline);
        vkCmdBindDescriptorSets(context.command_buffer[command_idx], VK_PIPELINE_BIND_POINT_COMPUTE,
                                pipeline.pipeline_layout, 0, descriptor_sets.size(), descriptor_sets.data(), 0, NULL);
        vkCmdPushConstants(context.command_buffer[command_idx], pipeline.pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT,
                           0, sizeof(UBOData), &data);
        const auto num_wg_x = (uint32_t)ceil(width / float(wg_x));
        const auto num_wg_y = (uint32_t)ceil(height / float(wg_y));
        vkCmdDispatch(context.command_buffer[command_idx], num_wg_x, num_wg_y, 1);
    }

    void create_fences() {
        VkFenceCreateInfo fence_create_info = {};
        fence_create_info.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
        fence_create_info.flags = 0;
        VK_CHECK(vkCreateFence(context.device, &fence_create_info, NULL, &fence));
    }

    void run_command_buffer(int command_idx) {

        VkSubmitInfo submit_info = {};
        submit_info.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
        submit_info.commandBufferCount = 1;
        submit_info.pCommandBuffers = &context.command_buffer[command_idx];
        VK_CHECK(vkQueueSubmit(context.compute_queue, 1, &submit_info, fence));
        VK_CHECK(vkWaitForFences(context.device, 1, &fence, VK_TRUE, 100000000000));
        vkResetFences(context.device, 1, &fence);
    }

    void cleanup(const std::vector<Pipeline> &pipelines) {
        if (enable_validation) {
            auto func =
                (PFN_vkDestroyDebugReportCallbackEXT)vkGetInstanceProcAddr(instance, "vkDestroyDebugReportCallbackEXT");
            if (func == nullptr) {
                throw std::runtime_error("Could not load vkDestroyDebugReportCallbackEXT");
            }
            func(instance, debug_report_callback, NULL);
        }

        for (auto pipeline : pipelines) {
            vkDestroyShaderModule(context.device, pipeline.shader_module, NULL);
        }
        vkDestroyDescriptorPool(context.device, descriptor_pool, NULL);
        for (int i = 0; i < MAX_DESCRIPTOR_SETS; i++) {
            vkDestroyDescriptorSetLayout(context.device, descriptor_set_layouts[i], NULL);
        }
        for (auto pipeline : pipelines) {
            vkDestroyPipelineLayout(context.device, pipeline.pipeline_layout, NULL);
        }
        for (auto pipeline : pipelines) {
            vkDestroyPipeline(context.device, pipeline.pipeline, NULL);
        }
        vkDestroyFence(context.device, fence, NULL);
        vkDestroyCommandPool(context.device, context.command_pool, NULL);
        vkDestroyDevice(context.device, NULL);
        vkDestroyInstance(instance, NULL);
    }
    Context context;
    UBOData data;

  private:
    VkInstance instance;
    VkDebugReportCallbackEXT debug_report_callback;
    // VkPipeline pipeline;
    // VkPipelineLayout pipeline_layout;
    // VkShaderModule compute_shader_module;
    VkDescriptorPool descriptor_pool;
    std::vector<VkDescriptorSet> descriptor_sets;
    std::vector<VkDescriptorSetLayout> descriptor_set_layouts;
    VkFence fence;
    // VkBuffer buffer;
    // VkDeviceMemory buffer_memory;

    // uint32_t buffer_size; // size of `buffer` in bytes.
    std::vector<const char *> enabled_layers;

    uint32_t compute_queue_family_index;
};

void barrier(GPUSimulation &simulation, Buffer &buffer, int command_idx,
             VkAccessFlags src_access = VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_SHADER_READ_BIT,
             VkAccessFlags dst_access = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT) {
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
    int num_wgs = ceil(dim / 1024.0f);
    vkCmdFillBuffer(simulation.context.command_buffer[command_idx], counter_buffer.handle, 0, counter_buffer.size, 1);
    fill_barrier = buffer_barrier(counter_buffer.handle, VK_ACCESS_TRANSFER_WRITE_BIT,
                                  VK_ACCESS_SHADER_WRITE_BIT | VK_ACCESS_SHADER_READ_BIT);
    vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_TRANSFER_BIT,
                         VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &fill_barrier, 0, 0);
    std::array<VkBufferMemoryBarrier, 2> barriers{res_barrier, fill_barrier};
    while (num_wgs != 1) {
        simulation.record_command_buffer(reduce_pipeline, command_idx, 1024, 1, num_wgs, 1);
        num_wgs = ceil(num_wgs / 1024.0f);
        if (num_wgs > 1) {
            vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                                 VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, barriers.data(), 0, 0);
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
    vkCmdPipelineBarrier(simulation.context.command_buffer[command_idx], VK_PIPELINE_STAGE_TRANSFER_BIT,
                         VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 1, &fill_barrier, 0, 0);
    std::array<VkBufferMemoryBarrier, 2> barriers{res_barrier, fill_barrier};
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

void scalar_div(GPUSimulation &simulation, Pipeline &pipeline, int command_idx) {
    simulation.record_command_buffer(pipeline, command_idx, 1, 1, 1, 1);
}

void vec_saxpy(GPUSimulation &simulation, Pipeline &pipeline, int command_idx, int dim) {
    simulation.record_command_buffer(pipeline, command_idx, 1024, 1, dim, 1);
}
