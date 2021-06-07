### Cases

a) [Plane shear flow](#a-plane-shear-flow)  
b) [The Karman Vortex Street](#b-the-karman-vortex-street)  
c) [Flow over a Step](#c-flow-over-a-step)  
d) [Natural Convection](#d-natural-convection)  
1. [nu = 0.001](#1-nu-0001)  
2. [nu = 0.0002](#2-nu-00002)  

e) [Fluid Trap](#e-fluid-trap)  
1. [T_h to T_c](#1-t_h-to-t_c)  
2. [T_c to T_h](#2-t_c-to-t_h)  

f) [Rayleigh-Bernard Convection](#f-rayleigh-bernard-convection)  

### a) Plane shear flow 
**Velocity**<br>
![PlaneShear/uv.png](./PlaneShear/uv.png)  
**Stream tracer**<br>
![PlaneShear/streamtracer.png](./PlaneShear/streamtracer.png)  
**Glyph**<br>
![PlaneShear/glyph.png](./PlaneShear/glyph.png)  
**Pressure**<br>
![PlaneShear/p.png](./PlaneShear/p.png)  
[Return to the cases](#cases)



## Performance Analysis  

### Strong Scaling on RayleighBenardConvection with

* imax = 85  
* jmax = 18  
* t_end = 20000
* 1-8 processes, domain only split along x-axis

**Speedup:**  
![StrongScaling/Speedup.png](./StrongScaling/Speedup.png)  
**Parallel Efficiency**  
![StrongScaling/Efficiency.png](./StrongScaling/Efficiency.png)  

For 2 processes the speedup can be really close to the theoretical speedup. There are multiple reasons for why the speedup stays far below the theoretical speedup for more processes. With increasing number of processes we increase also the communication overhead and increase the fraction of sequential code. The efficiency drops for more than 5 processes since then --use-hwthread-cpus needed to be enabled. The big drop for 8 processes is mainly due to our to simple partitioning in which the last processor takes additionally the leftover cells. Hence the work is not properly distributed for when the number of cells is not dividable by the number of processes.

### Weak Scaling on LidDrivenCavity with

* imax =     50 | 100 | 150 | 100 | 250 | 150 | 350 | 200
* jmax =     50 | 50  | 50  | 100 | 50  | 100 | 50  | 100  
* xlength =  1 | 2 | 3 | 2 | 5 | 3 | 7 | 8  
* ylength =  1 | 1 | 1 | 2 | 1 | 2 | 1 | 2

**Runtime:**  
![WeakScaling/Runtime.png](./WeakScaling/Runtime.png)  

For the weak scaling analysis we try to increase the domain with the number of processes to keep the workload per process approximately the same. Additionally to imax and jmax we adjusted xlength and ylength to avoid more timesteps due to smaller grid spacings.  