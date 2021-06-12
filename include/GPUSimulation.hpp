#include <vector>
#include <vulkan/vulkan.h>
#ifdef NDEBUG
constexpr bool enable_validation = true;
#else
constexpr bool enable_validation = true;
#endif

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

struct UBOData {
    float dt;
    float parity = 0;
    uint32_t fluid_cells_size;
};


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
    int idx = 0;
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
                VkSharingMode sharing_mode, VkDeviceSize size, void *data = nullptr, bool use_staging = false) {
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
                                  VK_SHARING_MODE_EXCLUSIVE, size, data);
            this->create(ctx, VK_BUFFER_USAGE_TRANSFER_DST_BIT | usage, mem_property_flags, sharing_mode, size, nullptr,
                         /* use_staging */ false);

            // @performance
            VK_CHECK(vkResetCommandPool(ctx->device, ctx->command_pool, 0));
            VkCommandBufferBeginInfo begin_info = {};
            begin_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
            begin_info.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
            VK_CHECK(vkBeginCommandBuffer(ctx->command_buffer[ctx->idx], &begin_info));
            VkBufferCopy copy_region = {0, 0, VkDeviceSize(size)};
            vkCmdCopyBuffer(ctx->command_buffer[ctx->idx], staging_buffer.handle, this->handle, 1, &copy_region);
            VkBufferMemoryBarrier copy_barrier =
                buffer_barrier(handle, VK_ACCESS_TRANSFER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT);
            vkCmdPipelineBarrier(ctx->command_buffer[ctx->idx], VK_PIPELINE_STAGE_TRANSFER_BIT,
                                 VK_PIPELINE_STAGE_ALL_COMMANDS_BIT, VK_DEPENDENCY_BY_REGION_BIT, 0, 0, 1,
                                 &copy_barrier, 0, 0);
            vkEndCommandBuffer(ctx->command_buffer[ctx->idx]);
            VkSubmitInfo submitInfo = {VK_STRUCTURE_TYPE_SUBMIT_INFO};
            submitInfo.commandBufferCount = 1;
            submitInfo.pCommandBuffers = &ctx->command_buffer[ctx->idx];

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
    void upload(VkDeviceSize size, void *data = nullptr, Buffer *staging_buffer = nullptr) {
        if ((mem_prop_flags & (VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT))) {
            memcpy(this->data, data, size);
        } else {
            // Assume device local memory
            staging_buffer->upload(size, data);
            VK_CHECK(vkResetCommandPool(ctx->device, ctx->command_pool, 0));
            VkCommandBufferBeginInfo begin_info = {};
            begin_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
            begin_info.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
            VK_CHECK(vkBeginCommandBuffer(ctx->command_buffer[ctx->idx], &begin_info));
            VkBufferCopy copy_region = {0, 0, VkDeviceSize(size)};
            vkCmdCopyBuffer(ctx->command_buffer[ctx->idx], staging_buffer->handle, this->handle, 1, &copy_region);
            VkBufferMemoryBarrier copy_barrier =
                buffer_barrier(handle, VK_ACCESS_TRANSFER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT);
            vkCmdPipelineBarrier(ctx->command_buffer[ctx->idx], VK_PIPELINE_STAGE_TRANSFER_BIT,
                                 VK_PIPELINE_STAGE_ALL_COMMANDS_BIT, VK_DEPENDENCY_BY_REGION_BIT, 0, 0, 1,
                                 &copy_barrier, 0, 0);
            vkEndCommandBuffer(ctx->command_buffer[ctx->idx]);
            VkSubmitInfo submitInfo = {VK_STRUCTURE_TYPE_SUBMIT_INFO};
            submitInfo.commandBufferCount = 1;
            submitInfo.pCommandBuffers = &ctx->command_buffer[ctx->idx];
            VK_CHECK(vkQueueSubmit(ctx->compute_queue, 1, &submitInfo, VK_NULL_HANDLE));
            VK_CHECK(vkDeviceWaitIdle(ctx->device));
        }
    }
    void copy(Buffer &dst_buffer, bool reset = true) {
        VkBufferCopy copy_region;
        copy_region.srcOffset = 0;
        copy_region.dstOffset = 0;
        copy_region.size = size;
        VkCommandBufferBeginInfo begin_info = {};
        begin_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
        begin_info.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
        if (reset) {
            VK_CHECK(vkResetCommandPool(ctx->device, ctx->command_pool, 0));
            VK_CHECK(vkBeginCommandBuffer(ctx->command_buffer[ctx->idx], &begin_info)); 
        }
      
        vkCmdCopyBuffer(ctx->command_buffer[ctx->idx], handle, dst_buffer.handle, 1, &copy_region);
        VkBufferMemoryBarrier copy_barrier =
            buffer_barrier(handle, VK_ACCESS_TRANSFER_WRITE_BIT, VK_ACCESS_SHADER_READ_BIT);
        vkCmdPipelineBarrier(ctx->command_buffer[ctx->idx], VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_ALL_COMMANDS_BIT,
                             VK_DEPENDENCY_BY_REGION_BIT, 0, 0, 1, &copy_barrier, 0, 0);
        if (reset) {
            vkEndCommandBuffer(ctx->command_buffer[ctx->idx]);
            VkSubmitInfo submitInfo = {VK_STRUCTURE_TYPE_SUBMIT_INFO};
            submitInfo.commandBufferCount = 1;
            submitInfo.pCommandBuffers = &ctx->command_buffer[ctx->idx];
            VK_CHECK(vkQueueSubmit(ctx->compute_queue, 1, &submitInfo, VK_NULL_HANDLE));
            VK_CHECK(vkDeviceWaitIdle(ctx->device)); 
        }
    }
};

struct Descriptor {
    Buffer handle;
    int binding;
    VkDescriptorType desc_type = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
};

class GPUSimulation {
  public:
    void init() {
        create_instance();
        find_physical_device();
        create_device();
        create_descriptor_pool();
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

    // Returns the index of a queue family that supports compute operations.
    uint32_t get_compute_queue_family_index() {

        uint32_t queue_family_count;
        vkGetPhysicalDeviceQueueFamilyProperties(context.physical_device, &queue_family_count, NULL);

        // Retrieve all queue families.
        std::vector<VkQueueFamilyProperties> queue_families(queue_family_count);
        vkGetPhysicalDeviceQueueFamilyProperties(context.physical_device, &queue_family_count, queue_families.data());

        // Now find a family that supports compute.
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
        // Create the logical device
        VkDeviceQueueCreateInfo queue_create_info = {};
        queue_create_info.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
        compute_queue_family_index = get_compute_queue_family_index(); // find queue family with compute capability.
        queue_create_info.queueFamilyIndex = compute_queue_family_index;
        queue_create_info.queueCount = 1; // create one queue in this family. We don't need more.
        float queue_priorities = 1.0;     // we only have one queue, so this is not that imporant.
        queue_create_info.pQueuePriorities = &queue_priorities;

        VkDeviceCreateInfo device_create_info = {};

        // None for now
        VkPhysicalDeviceFeatures device_features = {};

        device_create_info.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
        device_create_info.enabledLayerCount = enabled_layers.size();
        device_create_info.ppEnabledLayerNames = enabled_layers.data();
        device_create_info.pQueueCreateInfos =
            &queue_create_info; // when creating the logical device, we also specify what queues it has.
        device_create_info.queueCreateInfoCount = 1;
        device_create_info.pEnabledFeatures = &device_features;

        VK_CHECK(vkCreateDevice(context.physical_device, &device_create_info, NULL, &context.device));

        // Get a handle to the only member of the queue family.
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
        /*
              Our descriptor pool can only allocate a single storage buffer.
              */
        std::vector<VkDescriptorPoolSize> pool_sizes = {{VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 32},
                                                        {VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 32}};
        VkDescriptorPoolCreateInfo descriptor_pool_create_info = {};
        descriptor_pool_create_info.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
        descriptor_pool_create_info.maxSets = 32; // we only need to allocate one descriptor set from the pool.
        descriptor_pool_create_info.poolSizeCount = pool_sizes.size();
        descriptor_pool_create_info.pPoolSizes = pool_sizes.data();

        // create descriptor pool.
        VK_CHECK(vkCreateDescriptorPool(context.device, &descriptor_pool_create_info, NULL, &descriptor_pool));
    }

    void create_descriptor_set_layout(const std::vector<VkDescriptorSetLayoutBinding> &bindings) {
        VkDescriptorSetLayoutCreateInfo dsl_create_info = {VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
        dsl_create_info.bindingCount = bindings.size();
        dsl_create_info.pBindings = bindings.data();
        // Create the descriptor set layout.
        VK_CHECK(vkCreateDescriptorSetLayout(context.device, &dsl_create_info, NULL, &descriptor_set_layout));
        /*
      With the pool allocated, we can now allocate the descriptor set.
      */
        VkDescriptorSetAllocateInfo descriptor_set_allocate_info = {VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
        descriptor_set_allocate_info.descriptorPool = descriptor_pool; // pool to allocate from.
        descriptor_set_allocate_info.descriptorSetCount = 1;           // allocate a single descriptor set.
        descriptor_set_allocate_info.pSetLayouts = &descriptor_set_layout;

        // allocate descriptor set.
        VK_CHECK(vkAllocateDescriptorSets(context.device, &descriptor_set_allocate_info, &descriptor_set));
    }

    void push_descriptors(const std::vector<Descriptor> &descriptors) {
        // Specify the buffer to bind to the descriptor.
        std::vector<VkWriteDescriptorSet> write_descriptor_sets;
        for (const auto &descriptor : descriptors) {
            VkWriteDescriptorSet write_descriptor_set = {};
            write_descriptor_set.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            write_descriptor_set.dstSet = descriptor_set;         // write to this descriptor set.
            write_descriptor_set.dstBinding = descriptor.binding; // write to the first, and only binding.
            write_descriptor_set.descriptorCount = 1;             // update a single descriptor.
            write_descriptor_set.descriptorType = descriptor.desc_type;
            write_descriptor_set.pBufferInfo = &descriptor.handle.descriptor;
            write_descriptor_sets.push_back(write_descriptor_set);
        }
        // perform the update of the descriptor set.
        vkUpdateDescriptorSets(context.device, write_descriptor_sets.size(), write_descriptor_sets.data(), 0, NULL);

        /*  vkCmdBindDescriptorSets(command_buffer, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline_layout, 0, 1,
           &descriptor_set, 0, NULL);*/
    }

    // Read file into array of bytes, and cast to uint32_t*, then return.
    // The data has been padded, so that it fits into an array uint32_t.
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
        // Create shader module
        uint32_t filelength;
        // the code in comp.spv was created by running the command:
        // glslangValidator.exe -V shader.comp
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

        /*
        Now let us actually create the compute pipeline.
        A compute pipeline is very simple compared to a graphics pipeline.
        It only consists of a single stage with a compute shader.

        So first we specify the compute shader stage, and it's entry point(main).
        */

        VkPipelineShaderStageCreateInfo shader_stage_create_info = {};
        shader_stage_create_info.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        shader_stage_create_info.stage = VK_SHADER_STAGE_COMPUTE_BIT;
        shader_stage_create_info.module = compute_shader_module;
        if (specialization_data != -1) {
            shader_stage_create_info.pSpecializationInfo = &specialization_info;
        }
        shader_stage_create_info.pName = "main";

        /*
        The pipeline layout allows the pipeline to access descriptor sets.
        So we just specify the descriptor set layout we created earlier.
        */
        VkPipelineLayoutCreateInfo pipeline_layout_create_info = {};
        pipeline_layout_create_info.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
        pipeline_layout_create_info.setLayoutCount = 1;
        pipeline_layout_create_info.pSetLayouts = &descriptor_set_layout;
        VK_CHECK(vkCreatePipelineLayout(context.device, &pipeline_layout_create_info, NULL, &pipeline_layout));

        VkComputePipelineCreateInfo pipeline_create_info = {};
        pipeline_create_info.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
        pipeline_create_info.stage = shader_stage_create_info;
        pipeline_create_info.layout = pipeline_layout;
        VK_CHECK(vkCreateComputePipelines(context.device, VK_NULL_HANDLE, 1, &pipeline_create_info, NULL, &pipeline));
        return {pipeline, pipeline_layout, compute_shader_module};
    }

    void create_command_pool() {
        /*
              We are getting closer to the end. In order to send commands to the device(GPU),
              we must first record commands into a command buffer.
              To allocate a command buffer, we must first create a command pool. So let us do that.
              */
        VkCommandPoolCreateInfo command_pool_create_info = {};
        command_pool_create_info.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
        command_pool_create_info.flags = 0;
        // the queue family of this command pool. All command buffers allocated from this command pool,
        // must be submitted to queues of this family ONLY.
        command_pool_create_info.queueFamilyIndex = compute_queue_family_index;
        VK_CHECK(vkCreateCommandPool(context.device, &command_pool_create_info, NULL, &context.command_pool));
        VkCommandBufferAllocateInfo command_buffer_allocate_info = {};
        command_buffer_allocate_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
        command_buffer_allocate_info.commandPool = context.command_pool; // specify the command pool to allocate from.
        // if the command buffer is primary, it can be directly submitted to queues.
        // A secondary buffer has to be called from some primary command buffer, and cannot be directly
        // submitted to a queue. To keep things simple, we use a primary command buffer.
        command_buffer_allocate_info.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
        context.command_buffer.resize(3);
        command_buffer_allocate_info.commandBufferCount = 3; // allocate a single command buffer.
        VK_CHECK(vkAllocateCommandBuffers(context.device, &command_buffer_allocate_info,
                                          context.command_buffer.data())); // allocate command buffer.
    }

    void begin_recording(VkCommandBufferUsageFlags flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT) {
        /*
           Now we shall start recording commands into the newly allocated command buffer.
           */
        VK_CHECK(vkResetCommandPool(context.device, context.command_pool, 0));
        VkCommandBufferBeginInfo begin_info = {};
        begin_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
        begin_info.flags = flags; // the buffer is only submitted and used once
                                  // in
                                  // this application.
        VK_CHECK(vkBeginCommandBuffer(context.command_buffer[context.idx], &begin_info)); // start recording commands.
    }

    void end_recording() {
        VK_CHECK(vkEndCommandBuffer(context.command_buffer[context.idx])); // end recording commands.
    }

   void begin_end_record_command_buffer(Pipeline pipeline, int wg_x = 32, int wg_y = 32, int width = 102,
                                         int height = 22) {
        begin_recording();
        record_command_buffer(pipeline, wg_x, wg_y, width, height);
        end_recording();
    }

    void record_command_buffer(Pipeline pipeline, int wg_x = 32, int wg_y = 32, int width = 102, int height = 22) {
        /*
        We need to bind a pipeline, AND a descriptor set before we dispatch.

        The validation layer will NOT give warnings if you forget these, so be very careful not to forget them.
        */
        vkCmdBindPipeline(context.command_buffer[context.idx], VK_PIPELINE_BIND_POINT_COMPUTE, pipeline.pipeline);

        vkCmdBindDescriptorSets(context.command_buffer[context.idx], VK_PIPELINE_BIND_POINT_COMPUTE, pipeline.pipeline_layout, 0, 1,
                                &descriptor_set, 0, NULL);

        const auto num_wg_x = (uint32_t)ceil(width / float(wg_x));
        const auto num_wg_y = (uint32_t)ceil(height / float(wg_y));
        vkCmdDispatch(context.command_buffer[context.idx], num_wg_x, num_wg_y, 1);
    }

    void create_fences() {
        VkFenceCreateInfo fence_create_info = {};
        fence_create_info.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
        fence_create_info.flags = 0;
        VK_CHECK(vkCreateFence(context.device, &fence_create_info, NULL, &fence));
    }

    void run_command_buffer() {

        VkSubmitInfo submit_info = {};
        submit_info.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
        submit_info.commandBufferCount = 1;
        submit_info.pCommandBuffers = &context.command_buffer[context.idx];
        VK_CHECK(vkQueueSubmit(context.compute_queue, 1, &submit_info, fence));
        /*
        The command will not have finished executing until the fence is signalled.
        So we wait here.
        We will directly after this read our buffer from the GPU,
        and we will not be sure that the command has finished executing unless we wait for the fence.
        Hence, we use a fence here.
        */
        VK_CHECK(vkWaitForFences(context.device, 1, &fence, VK_TRUE, 100000000000));
        vkResetFences(context.device, 1, &fence);
    }

    void cleanup(const std::vector<Pipeline> &pipelines) {
        /*
        Clean up all Vulkan Resources.
        */
        if (enable_validation) {
            // destroy callback.
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
        vkDestroyDescriptorSetLayout(context.device, descriptor_set_layout, NULL);
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

  private:
    VkInstance instance;
    VkDebugReportCallbackEXT debug_report_callback;
    // VkPipeline pipeline;
    // VkPipelineLayout pipeline_layout;
    // VkShaderModule compute_shader_module;
    VkDescriptorPool descriptor_pool;
    VkDescriptorSet descriptor_set;
    VkDescriptorSetLayout descriptor_set_layout;
    VkFence fence;
    // VkBuffer buffer;
    // VkDeviceMemory buffer_memory;

    // uint32_t buffer_size; // size of `buffer` in bytes.
    std::vector<const char *> enabled_layers;

    uint32_t compute_queue_family_index;
};