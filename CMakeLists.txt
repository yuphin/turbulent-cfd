# Definition of the minimum required cmake Version
cmake_minimum_required(VERSION 3.0)
# Definition of the Project
# Later you can access the project variable like ${CFDLAB_SOURCE_DIR}
project(CFDLAB VERSION 1.0 LANGUAGES CXX)

option(USE_VULKAN "USE_VULKAN" OFF)
option(USE_CUDA "USE_CUDA" OFF)

file(GLOB_RECURSE SRCS src/*.cpp)
set(CUDA_SRC src/CUDASolver.cu)

# Definition of the C++ Standard 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# Require a package
find_package(MPI REQUIRED)
# find_package(Boost COMPONENTS filesystem REQUIRED)

# VTK Library
if(WIN32)
  find_package(VTK REQUIRED)
else()
  find_package(VTK 7 REQUIRED)
endif()
if (NOT VTK_FOUND)
  message("VTK not found")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (${VTK_VERSION} VERSION_LESS 9)
  include(${VTK_USE_FILE})
endif()

# Filesystem
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 9)
    message("g++ Version is 9 or higher")
    add_compile_definitions(gpp9)
else()
    message("g++ Version is lower than Version 9")
    if (NOT(APPLE OR WIN32))  
      # Link Filesystem Libraries
      link_libraries(stdc++fs)
    endif()
endif()

# Add include directory
include_directories(include)

if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX /usr/local )
endif()

if(USE_VULKAN)
    find_package(Vulkan)
    add_definitions(-DUSE_VULKAN)
    include_directories(${Vulkan_INCLUDE_DIR})
else()
    list(REMOVE_ITEM SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/VulkanUtils.cpp)
    list(REMOVE_ITEM SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/VulkanSolver.cpp)
endif()

if(USE_CUDA)
    set(CMAKE_CUDA_ARCHITECTURES 86)
    enable_language(CUDA)
    add_executable(vortigen ${SRCS} ${CUDA_SRC})
    add_definitions(-DUSE_CUDA)
    set_target_properties(vortigen PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    include_directories(${CUDA_INCLUDE_DIRS})
    include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
else()
     add_executable(vortigen ${SRCS})

endif()

target_link_libraries(vortigen PRIVATE MPI::MPI_CXX)
target_link_libraries(vortigen PRIVATE ${VTK_LIBRARIES})
if(USE_VULKAN)
target_link_libraries(vortigen PRIVATE ${Vulkan_LIBRARY})
endif()
if(USE_CUDA)
target_link_libraries(vortigen ${CUDA_RUNTIME_LIBRARY})
endif()

install(TARGETS vortigen DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
