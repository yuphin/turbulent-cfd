# Definition of the minimum required cmake Version
cmake_minimum_required(VERSION 3.0)
# Definition of the Project
# Later you can access the project variable like ${CFDLAB_SOURCE_DIR}
project(CFDLAB VERSION 1.0)

# Definition of the C++ Standard 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_definitions(-Dsolution_liddriven)
add_definitions(-Dsolution_energy)
add_definitions(-Dsolution_parallelization)

# You can find package likes
# find_package(MPI)
# Require a package
find_package(MPI REQUIRED)
# Find a package with different components e.g. BOOST
# find_package(Boost COMPONENTS filesystem REQUIRED)

# VTK Library
find_package(VTK REQUIRED)
if (NOT VTK_FOUND)
  message("VTK not found")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
include(${VTK_USE_FILE})

# Filesystem
set(gpp9 False)
#set(gpp9 True)
if(gpp9)
    message("g++ Version is 9 or higher")
    add_compile_definitions(gpp9)
else()
    message("g++ Version is lower than Version 9")
    if (NOT APPLE)  
      # Link Filesystem Libraries
      link_libraries(stdc++fs)
    endif()
endif()

# Add include directory
include_directories(include)

if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX /usr/local )
endif()

# Creating the executable of our project and the required dependencies
# the executable is called fluidchen
add_executable(fluidchen src/main.cpp src/Boundary.cpp src/Case.cpp src/Cell.cpp src/Discretization.cpp src/Fields.cpp src/Grid.cpp src/PressureSolver.cpp src/Communication.cpp)

# if you use external libraries you have to link them like 
target_link_libraries(fluidchen PRIVATE MPI::MPI_CXX)
target_link_libraries(fluidchen PRIVATE ${VTK_LIBRARIES})

install(TARGETS fluidchen DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
# If you write tests, you can include your subdirectory (in this case tests) as done here
# Testing
