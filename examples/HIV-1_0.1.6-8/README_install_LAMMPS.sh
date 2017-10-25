# Step 0: Prepare your system for LAMMPS
#
# LAMMPS requires MPI and FFTW
#        (Additional GPU-specific libraries are optional)
# If you are using ubuntu or windows WSL, install the debian packages you need:
# sudo apt-get install build-essential mpi-default-bin mpi-default-dev
# sudo apt-get install libfftw3-dev libjpeg-dev libpng12-dev


# STEP 1: Download LAMMPS
git clone -b master https://github.com/lammps/lammps.git lammps

# STEP 2: Select what features of LAMMPS you need.  I recommend the following:
cd lammps/src
make yes-molecule
make yes-rigid
make yes-kspace
make yes-manybody
make yes-mc
make yes-replica
make yes-body
make yes-asphere
make yes-class2
make yes-dipole
make yes-user-misc
#make yes-kokkos  # <- GPU acceleration  (avoid "yes-cuda" or "yes-gpu")

# STEP 3: Compile LAMMPS

make ubuntu

# This should eventually create a file named "lmp_ubuntu".
# 
# STEP 4: Copy this file to somewhere in your PATH
#
# Alternatively, if you use a different flavor of linux, install mpi and fftw3
# and use:
#
# make mpi
#
# This should create a file named "lmp_mpi".
