# 3D-flat-bystander

Three-species bystander model

# simulation.cc:
Performs a single run simulation of the population growth and prints the coordinates of cells (x,y,z) and the strain of the cells (1, 2, or 3) to an output file named "simulation_runXX"

The code is run with input parameters in the command line as: "./simulation LATTSIZE NGEN"
Where LATTSIZE is the size of one side of the square lattice (in the spatial direction) and NGEN is the number of generations (the temporal direction).

# phase.cc:
Calculates the density of strains at the population front after some number of generations (NGEN) for varying parameter values {mu,s} (mutation and selection). Prints the densities for each pair of values to an output file named "phase_runXX"

The code is run with input parameters in the command line as: "./simulation LATTSIZE NGEN NSLICES NRUNS"
Where LATTSIZE and NGEN is as before. NSLICES is the number of sampled points in the parameter space, and NRUNS is the number of runs to average over. This code can be implemented on any number of computing cores using mpirun.
