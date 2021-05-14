# Design of boundaries and cell types

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

### Information for implementing Boundaries

In accordance with the cell types there are 4 boundary classes.  
While iterating over the cells use the cells id to get information on velocity and temperature of the boundary.

```C++
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        int id = cell->id();
        
        // For walls
        double velocity = _wall_velocity[Id];
        double temperature = _wall_temperature[Id];
        // For inlets
        double u_velocity = _inlet_U[Id];
        double v_velocity = _inlet_V[Id];
        double v_velocity = _inlet_V[Id];
```
