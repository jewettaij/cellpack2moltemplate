# The default name of the LAMMPS binary depends on the operating system
# where you installed it.
# These instructions assume you have installed LAMMPS on an ubuntu linux
# distribution (or something related such as mint or debian, or windows WSL).
# If you are using a different flavor of linux, then you likely will want to
# replace the "lmp_ubuntu" commands with "lmp_mpi" 



# There are 2 steps required for running most LAMMPS simulations:


# 1) Minimize the system.
# Regardless how careful you are, some of the atoms will overlap slightly.
# Before running a simulation, you must relax the system.
# There are multiple ways to do this.  
# One strategy outlined in the "run.in.min1" and "run.in.min2" files


# There are two steps in the process of minmization in this example.
# LAMMPS must be run separately for each step:
lmp_ubuntu -i run.in.min1  # <-- use soft "gaussian" forces + annealing
lmp_ubuntu -i run.in.min2  # <-- use hard "lennard jones" forces + langevin

# or alternately, to run LAMMPS in parallel, try:
# mpiexec -np 4 lmp_ubuntu -i run.in.min1
# mpiexec -np 4 lmp_ubuntu -i run.in.min2

# (Minimizing the HIV example with ~100000 particles takes approximately 
#  45 minutes or so on a 2-core laptop without GPU acceleration.)
# This will create a file named "system_after_min.data" which we will
# need in the next step.
# Note during the minimization process, some of the proteins fly out of the core
# This can be fixed by improving the minimization protocol or tinkering
# with the lammps "minimize" code.  (Let's worry about this later.)





# 2) Run the main simulation


lmp_ubuntu -i run.in

# or
# mpiexec -np 4 lmp_ubuntu -i run.in


#    (Currently the system disolves since there are no forces to
#     hold things together.)
