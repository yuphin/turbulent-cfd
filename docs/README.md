### Introduction
In this CFD project, our main goal was to implement and investigate different turbulent models. Aside from this, we wanted to explore different solvers, parallelization of our code on the GPU with different APIs, preconditioning, adaptive grid refinement and more. Our CFD solver comprises of 3 types of simulations: a CPU, a Vulkan and a CUDA simulation. It also supports MPI on the CPU path. During the course of our project we have implemented K-epsilon, K-omega and K-omega SST turbulent models.

### Running
The code has been tested on both Linux and Windows. For Linux you can build the code via CMake.
Make sure you enable CUDA or Vulkan support from the CMake settings.
```
Vulkan : cmake .. -DUSE\_VULKAN=ON
CUDA : cmake .. -DUSE\_CUDA=ON
```
Otherwise simply run with
```
mkdir build
cd build
cmake ..
make
```
Make sure you delete you CMake cache if you change these settings midway through the build. 
Also note that the Vulkan code is only tested and designed on a RTX 3060 GPU, therefore there may be problems on lower-end hardware.

In order to run **Fluidchen**, the case file should be given as input parameter. Some default case files are located in the `example_cases` directory. If you installed **Fluidchen**, you can execute them from anywhere you want as  
* For Serial

```shell
fluidchen /path/to/case/case_name.dat [-log]
```
* For MPI
```
mpirun -np <num_processes> ./fluidchen /path/to/case/case_name.dat [-log] 
```

This will run the case file and create the output folder `/path/to/case/case_name_Output` which holds the `.vtk` files of the solution. If the `-log` flag is specified a log file will also be created in the output directory. The output folder is created in the same location as your case file. Note that this may require write permissions in the given directory.

If input file does not contain a geometry file, fluidchen will run lid-driven cavity case with given parameters.

### Scene configuration
In addition to the usual options from the previous exercises, we have other options. All properties relevant for the simulation can be set in a .dat file.  
For reference see [this template](TEMPLATE.dat)

* model
  ```
  0 -> Turbulence(default)
  1 -> K-epsilon model
  2 -> K-omega model
  3 -> K-omega SST model  
  ```
* solver
  ```
  0 -> SOR(default)
  1 -> PCG
  ```
* simulation
  ```
  0 -> CPU(default)
  1 -> Cuda
  2 -> Vulkan
  ```
* preconditioner
  ```
  -1 -> off(default)
   0 -> AINV
   1 -> SSOR
   2 -> Jacobi preconditioner
   ```
Note that following combinations are not supported:
  * PCG - MPI
  * Cuda - MPI
  * Vulkan - MPI

### Design of boundaries and cell types

There are 5 different types of cells:

* FLUID
* OUTLET
* INLET
* NOSLIP_WALL
* FREESLIP_WALL

Cell types and geometry can be set via .pgm files, where the following values correspond to respective cell types.

| Values    | Cell Type     |
| --------- |:-------------:|
| 0         | FLUID         |
| 1         | OUTLET        |
| 2 - 9     | INLET         |
| 10 - 19   | NOSLIP_WALL   |
| 20 - 25   | FREESLIP_WALL |


### Single vs Double Precision

By default, we use double precision floating point numbers. If you want to switch to single precision, use the macro 
```C++
   #define USE_FLOATS 1
```
in `include\Utilities.hpp`
