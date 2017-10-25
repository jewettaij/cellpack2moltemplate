# Step 0: Make sure you have installed moltemplate and cellpack2moltemplate
#
# git clone https://github.com/jewettaij/moltemplate moltemplate 
# git clone https://github.com/jewettaij/cellpack2moltemplate cellpack2moltemplate 
# # optional: use virtualenv to avoid polluting your python package system
# virtualenv venv_moltemplate
# source venv_moltemplate/bin/activate
#
# cd moltemplate
# pip install .
# cd ..
# cd cellpack2moltemplate
# pip install .
# cd ..


# Step 1: Convert the JSON file created by CellPACK into moltemplate format:

cellpack2lt.py -in HIV-1_0.1.6-8_mixed_radii_pdb.cpr -out system.lt

# (NOTE: The HIV-1_0.1.6-8_mixed_radii_pdb.cpr is too big to include here)

# Step 2: Convert the "system.lt" file into files that LAMMPS can read:
#         (This takes about 5 minutes on a ~100000 particle system.)

moltemplate.sh system.lt

# If you have VMD installed, then I recommend using this command instead:
# moltemplate.sh -vmd system.lt
# It will display the system and create a "system.psf" file you can use
# later to display the simulation trajectories you create with LAMMPS.
