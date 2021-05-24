## Running

In order to run **Fluidchen**, the case file should be given as input parameter. Some default case files are located in the `example_cases` directory. If you installed **Fluidchen**, you can execute them from anywhere you want as
For Serial:

```shell
fluidchen /path/to/case/case_name.dat [-log]
```

This will run the case file and create the output folder `/path/to/case/case_name_Output` which holds the `.vtk` files of the solution. If the `-log` flag is specified a log file will also be created in the output directory. The output folder is created in the same location as your case file. Note that this may require write permissions in the given directory.

If input file does not contain a geometry file, fluidchen will run lid-driven cavity case with given parameters.

## Design of boundaries and cell types

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

All properties relevant for the simulation can be set in a .dat file.  
For reference see [this template](TEMPLATE.dat)


## Single vs Double Precision

By default, we use double precision floating point numbers. If you want to switch to single precision, use the macro 
```C++
   #define USE_FLOATS 1
```
in `include\Utilities.hpp`
