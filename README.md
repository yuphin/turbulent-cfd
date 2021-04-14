![](FluidchenLogo.png)

Fluidchen is a CFD Solver developed for the CFD Lab taught at TUM Informatics, Chair of Scientific Computing in Computer Science.

Lecturers:
- Gerasimos Chourdakis, M.Sc.
- Benjamin Rüth, M.Sc.

Contributors:
- Oguz Ziya Koseomur
- Carme Homs Pons
- Aristotelis Economides
- Moritz Gnisia

We thank https://instagram.com/radioactivebunnies for designing such a cool logo for us:)

## Software Requirements

* VTK 7 or higher
* GCC 9 (optional) 
  
Detailed information is given below.

## Installing
```shell
git clone git@gitlab.lrz.de:YOUR_ACCOUNT_NAME/YOUR_GROUP_NAME.git
cd YOUR_GROUP_NAME
mkdir build && cd build
cmake ..
make
make install
```
These commands will create the executable `fluidchen` and copy it to the default directory `/usr/local/bin` . If you want to install to a custom path, execute the cmake command as
```shell
cmake -DCMAKE_INSTALL_PREFIX=/path/to/directory ..
```
After `make && make install` **fluidchen** will be installed to `/path/to/directory/bin` . Note that you may need to update your `PATH` environment variable.

By default, **fluidchen** is installed in `DEBUG` mode. To obtain full performance, you can execute cmake as
```shell
cmake -DCMAKE_BUILD_TYPE=RELEASE ..
```
or
```shell
cmake -DCMAKE_CXX_FLAGS="-O3" ..
```

## Running
In order to run **Fluidchen**, the case file should be given as input parameter. Some default case files are located in the `example_cases` directory. If you installed **Fluidchen**, you can execute them from anywhere you want as
For Serial:
```shell
fluidchen /path/to/fluidchen/example_cases/LidDrivenCavity/LidDrivenCavity.dat
```
This will run the case file and create the output folder `/path/to/case/case_name_Output` which holds the `.vtk` files of the solution. The output folder is created in the same location as your case file. Note that this may require write permissions in the given directory.

If input file does not contain a geometry file, fluidchen will run lid-driven cavity case with given parameters.

### GCC version

You can get you current version of GCC by running:

```shell
g++ -v
```

### Defining your GCC version

If you have GCC 9 or newer, you can set in the `CMakeLists.txt` file:

```cmake
set(gpp9 True)
```

If you have a version lower than 9, then you don't have to modify the `CMakeLists.txt` file.

This will affect how we are using the C++ filesystem library, which is available already in GCC 7 as an experimental feature.

### Setup of VTK and GCC 9 (Ubuntu **20.04**)

```
apt-get update &&
apt-get upgrade -y &&
apt-get install -y build-essential cmake libvtk7-dev libfmt-dev
```

### Setup of VTK and GCC 9 (Ubuntu **18.04**)

If you want, you can upgrade your compiler version to have access to more recent C++ features.
This is, however, optional.

```
apt-get update &&
apt-get install -y software-properties-common &&
add-apt-repository -y ppa:ubuntu-toolchain-r/test &&
apt-get upgrade -y &&
apt-get install -y build-essential cmake libvtk7-dev libfmt-dev gcc-9 g++-9
apt-get install -y gcc-9 g++-9
```

### Setup of VTK (OSX)
In macOS you do not need to install GCC compiler, you can directly use Clang (no need to do anything). The dependencies can be installed by using `homebrew` package manager as
```
brew install vtk
```

If you want to install GCC compiler anyways
````
brew install vtk
````

By default, `g++` command is linked to `clang` command in OSX operating systems. In order to use GCC
```shell
cmake -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-9 ..
```
should be used.

## Using CMake

CMake is a C++ build system generator, which simplifies the building process compared e.g. to a system-specific Makefile. The CMake configuration is defined in the `CMakeList.txt` file.

In order to build your code with CMake, you can follow this (quite common) procedure:

1. Create a build directory: `mkdir build`
2. Get inside it: `cd build`
3. Configure and generate the build system: `cmake ..` (Note the two dots, this means that the `CmakeLists.txt` File is in the folder above)
4. Build your code: `make` (build the executable)

### Troubleshooting: VTK not found

You might run into a problem where the VTK library is not found. To fix this, you can try the following steps:

1. Find the installation path of your VTK library 
2. Define this path as an environment variable, as e.g. `export VTK_DIR=".../lib/cmake/vtk-8.2"`
3. Start in a clean build folder
4. Run `cmake ..` again

### Set a different GCC version

If you have multiple compiler versions installed you can set the GCC version which should be used by `cmake` like this:

```shell
export CXX=`which g++-7`
```

Make sure to use a backtick (\`) to get the `which` command executed. Afterwards, you can run `cmake ..`.

