#!/usr/bin/env python

# Author: Andrew Jewett (jewett.aij at g mail)
# License: 3-clause BSD License  (See LICENSE.TXT)
# Copyright (c) 2017, California Institute of Technology
# All rights reserved.

"""
cellpack2lt.py converts json formatted files created by CellPACK into
moltemplate format.
"""

g_program_name = __file__.split('/')[-1]   # = 'cellpack2lt.py'
__version__ = '0.0.4'
__date__ = '2017-10-26'

g_control_vmd_colors = False


doc_msg = \
    'Typical Usage:\n\n' + \
    '  ' + g_program_name + ' -in HIV-1_0.1.cpr -out system.lt\n' + \
    ' or \n' + \
    '  ' + g_program_name + '  <  HIV-1_0.1.cpr  >   system.lt\n' + \
    '\n' + \
    'where \"HIV-1_0.1.cpr\" is a JSON file in CellPACK output format,\n' + \
    '  and \"HIV-1_0.1.lt\" is the corresponding file converted to moltemplate format\n' + \
    'Optional Arguments\n' + \
    '   -in FILE_NAME     # Specify the name of the input file (as opposed to \"<\")\n' + \
    '   -url URL          # Read the CellPACK JSON text from URL instead of a file\n' +\
    '   -out FILE_NAME    # Specify the output file (as opposed to \">\")\n' + \
    '   -pairstyle STYLE  # Select force formula (eg. -pairstyle lj/cut/coul/debye)\n' + \
    '   -debye            # Specify the Debye length (if applicable)\n' +\
    '   -epsilon          # Specify the \"Epsilon\" Lennard-Jones coeff (default: 1)\n' + \
    '   -deltaR           # Specify the resolution of particle radii (default 0.1)\n' + \
    '   -name OBJECTNAME  # Create a moltemplate object which contains everything.\n' + \
    '                     # (useful if you want multiple copies of the system later)\n'

import sys, json
from collections import defaultdict
import numpy as np


class InputError(Exception):
    """ A generic exception object containing a string for error reporting.
        (Raising this exception implies that the caller has provided
         a faulty input file or argument.)

    """

    def __init__(self, err_msg):
        self.err_msg = err_msg

    def __str__(self):
        return self.err_msg

    def __repr__(self):
        return str(self)


def Quaternion2Matrix(q, M):
    "convert a quaternion q to a 3x3 rotation matrix M"""

    M[0][0] =  (q[0]*q[0])-(q[1]*q[1])-(q[2]*q[2])+(q[3]*q[3])
    M[1][1] = -(q[0]*q[0])+(q[1]*q[1])-(q[2]*q[2])+(q[3]*q[3])
    M[2][2] = -(q[0]*q[0])-(q[1]*q[1])+(q[2]*q[2])+(q[3]*q[3])
    M[0][1] = 2*(q[0]*q[1] - q[2]*q[3]);
    M[1][0] = 2*(q[0]*q[1] + q[2]*q[3]);
    M[1][2] = 2*(q[1]*q[2] - q[0]*q[3]);
    M[2][1] = 2*(q[1]*q[2] + q[0]*q[3]);
    M[0][2] = 2*(q[0]*q[2] + q[1]*q[3]);
    M[2][0] = 2*(q[0]*q[2] - q[1]*q[3]);



def AffineTransformQ(x_new, x, q, deltaX, M=None):
    """ This function performs an affine transformation on vector "x".
        Multiply 3-dimensional vector "x" by first three columns of 3x3 rotation
        matrix M.  Add to this the final column of M.  Store result in "dest":
    x_new[0] = M[0][0]*x[0] + M[0][1]*x[1] + M[0][2]*x[2]  +  deltaX[0]
    x_new[1] = M[1][0]*x[0] + M[1][1]*x[1] + M[1][2]*x[2]  +  deltaX[1]
    x_new[2] = M[2][0]*x[0] + M[2][1]*x[1] + M[2][2]*x[2]  +  deltaX[2]
        where the rotation matrix, M, is determined by the quaternion, q.
        Optionally, the user can supply a preallocated 3x3 array (M)

    """
    if M == None:
        M = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    assert((M.shape[0] == 3) and (M.shape[1] == 3))

    Quaternion2Matrix(q, M)

    # There's probably a faster way to do this in numpy without using a for loop
    # but I'm too lazy to lookup what it is
    for i in range(0, 3):
        x_new[i] = 0.0
        for j in range(0, 3):
            x_new[i] += M[i][j] * x[j]
        x_new[i] += deltaX[i]  # (translation offset stored in final column)



def AdjustBounds(bounds, X):
    """ 
    Check if a coordinate lies within the box boundary.
    If not, update the box boundary.
    """
    assert(len(bounds) == 2)
    assert(len(bounds[0]) == len(bounds[1]) == len(X))

    for d in range(0, len(X)):
        if bounds[1][d] < bounds[0][d]:
            bounds[0][d] = bounds[1][d] = X[d]
        else:
            if X[d] < bounds[0][d]:
                bounds[0][d] = X[d]
            if bounds[1][d] < X[d]:
                bounds[1][d] = X[d]



def RadiiNeeded(tree,
                ir_needed,
                delta_r):
    """
    Note: LAMMPS (and most other) molecular dynamics simulation programs does
          not allow users to define unique radii for each of the 10^6 or so
          atoms in a simulation.  LAMMPS only allows for about 10^3 different
          -types- of particles in a simulation.  Each of the particle types
          in a CellPack file potentially has a unique radius, so we quantize
          each radius (by dividing it by delta_r and rounding) so that
          particles with similar radii are assigned to the same type.
          The first step is to figure out how many different radii of particles
          will be needed.  Do that by scanning the JSON file.

    This function recursively searches for objects in the JSON file containing
    a 'radii' field.  Quantize each of the radii (by dividing by delta_r) and 
    add them to the set of radii that we will need to define later (ir_needed).
    """

    if not isinstance(tree, dict):
        return
    if 'radii' in tree:
        r_ni = tree['radii']
        for n in range(0, len(tree['radii'])): #loop over "subunits" of this molecule
            for i in range(0, len(r_ni[n])):   #loop over atoms
                iradius = int(round(r_ni[n][i]/delta_r))  #(quantize the radii)
                ir_needed.add(iradius)
    else:
        for object_name in tree:
            RadiiNeeded(tree[object_name], ir_needed, delta_r)




def ConvertMolecule(molecule,
                    name,
                    l_mol_defs,
                    l_instances,
                    delta_r,
                    bounds):
    """ 
    Convert all inforation concering a type of molecule defined in a 
    CellPACK JSON file into moltemplate format.
    In this context, a \"molecule\" is defined as the smallest possible 
    subunit in the JSON files created by CellPACK, and it is assumed to
    contain 'positions', 'radii', and 'results' fields.
    Each molecule_node contains a definition of the type of molecule
    (these details are in the 'name', 'positions', 'radii' fields),
    as well as a list of separate instances (copies) of that molecule
    (in its 'results' field, which specifies the position and orientation
    of every copy of that type of molecule)
    This function converts this information into a moltemplate molecule
    objects, and a list of \"new\" commands which make copies of these
    molecules and move them into the correct positions and orientations.
    """

    l_mol_defs.append(name + ' inherits .../ForceField {\n')
    crds = []
    radii = []
    if not ('positions' in molecule):
        raise InputError('Error: Missing \"positions\" field in \"'+name+'\" molecule.\n')
    if not ('radii' in molecule):
        raise InputError('Error: Missing \"radii\" field in \"'+name+'\" molecule.\n')
    # Typically each molecule's "positions" and "radii" have an additional
    # level of hierachy beneath them. I don't know how many of them there are.
    # Loop over all of these "subunits" and collect all the coordinates
    # from each of them, appending the coordinates together into one big list:
    Xnid = molecule['positions']
    for n in range(0, len(molecule['positions'])):  #loop over "subunits" of this molecule
        for i in range(0, len(Xnid[n])):  #loop over atoms
            crds.append((Xnid[n][i][0],
                         Xnid[n][i][1],
                         Xnid[n][i][2]))
    # Do the same for the radii
    r_ni = molecule['radii']
    for n in range(0, len(molecule['radii'])): #loop over "subunits" of this molecule
        for i in range(0, len(r_ni[n])):       #loop over atoms
            radii.append(r_ni[n][i])
            #iradius = int(round(radii[i]/delta_r))  #(quantize the radii)
            #ir_needed.add(iradius)
    if len(crds) != len(radii):
        raise InputError('Error: \"positions\" and \"radii\" arrays in \"'+name+'\" have unequal length.\n')
    l_mol_defs.append('\n'
                      '      # AtomID MoleculeID AtomType Charge   X    Y    Z\n'
                      '\n')
    l_mol_defs.append('  write(\"Data Atoms\") {\n')
    assert(len(crds) == len(radii))
    for i in range(0, len(crds)):
        iradius = int(round(radii[i]/delta_r))  #(quantize the radii)
        atype_name = '@atom:A' + str(iradius)    #atom type depends on radius
        charge = '0.0'
        l_mol_defs.append('    $atom:a'+str(i+1)+'  $mol  '+atype_name+'  '+charge+'  '+str(crds[i][0])+' '+str(crds[i][1])+' '+str(crds[i][2])+'\n')
    l_mol_defs.append('  }  # end of: write(\"Data Atoms\") {...\n\n')
    l_mol_defs.append('}  # end of: \"'+name+'\" molecule definition\n\n')



    if not ('results' in molecule):
        raise InputError('Error: Missing \results\" field in \"'+name+'\" molecule.\n')

    deltaXs = []
    quaternions = []
    instances = molecule['results']
    for i in range(0, len(instances)):
        #if not (('0' in instance[i]) and ('1' in instance[i])):
        if not len(instances[i]) == 2:
            raise InputError('Error: Incorrect format in \results\" section of \"'+name+'\" molecule.\n')
        deltaXs.append((instances[i][0][0],
                        instances[i][0][1],
                        instances[i][0][2]))
        quaternions.append((instances[i][1][0],
                            instances[i][1][1],
                            instances[i][1][2],
                            instances[i][1][3]))

    l_instances.append('\n\n')
    #l_instances.append('# List of \"'+name+'\" instances:\n')
    for i in range(0, len(instances)):
        l_instances.append(name + '_instances[' + str(i) + '] = new ' + name +
                           '.quat(' + 
                           str(quaternions[i][0]) + ',' +
                           str(-quaternions[i][1]) + ',' +
                           str(-quaternions[i][2]) + ',' +
                           str(-quaternions[i][3]) +
                           ').move(' + 
                           str(deltaXs[i][0]) + ',' +
                           str(deltaXs[i][1]) + ',' +
                           str(deltaXs[i][2]) + ')\n')

        # Now determine the minimum/maximum coordinates of this object
        # and adjust the simulation box boundaries if necessary
        Xnid = molecule['positions']
        for n in range(0, len(molecule['positions'])):
            X_orig = [0.0, 0.0, 0.0]
            X = [0.0, 0.0, 0.0]
            for I in range(0, len(Xnid[n])):  #loop over atoms
                for d in range(0, 3):
                    X_orig[d] = molecule['positions'][n][I][d]
                    AffineTransformQ(X, X_orig,
                                     [quaternions[i][0],
                                      -quaternions[i][1],
                                      -quaternions[i][2],
                                      -quaternions[i][3]],
                                      deltaXs[i])
                AdjustBounds(bounds, X)

    if g_control_vmd_colors:
        for i in range(0, len(instances)):
            l_instances.append('\n')
            l_instances.append('write("vmd_commands.tcl") {  #(optional VMD file)\n')
            for I in range(0, len(radii)):
                r = iradius * delta_r
                atomid_name = '${atom:'+name+'_instances['+str(i)+']/a'+str(I+1)+'}'
                color_name  = '@color:'+name
                l_instances.append('  set sel [atomselect top "index '+atomid_name+'"]\n')
                l_instances.append('  \$sel set name '+color_name+'\n')
            l_instances.append('}  # end of "vmd_commands.tcl"\n')
            l_instances.append('\n')



def ConvertMolecules(molecules,
                     file_out,
                     delta_r,
                     bounds,
                     nindent=0):
    l_mol_defs = []
    l_instances = []
    for molecule_type_name in molecules:
        ConvertMolecule(molecules[molecule_type_name],
                        molecule_type_name,
                        l_mol_defs,
                        l_instances,
                        delta_r,
                        bounds)
        #when debugging, uncomment the next line:
        #break

    file_out.write('\n' + 
                   (nindent*'  ') + '# ----------- molecule definitions -----------\n'
                   '\n')
                   
    file_out.write((nindent*'  ') + (nindent*'  ').join(l_mol_defs))
    file_out.write('\n' +
                   nindent*'  ' + '# ----------- molecule instances -----------\n' +
                   '\n')
    file_out.write((nindent*'  ') + (nindent*'  ').join(l_instances))
    file_out.write('\n')




def ConvertSystem(tree,
                  file_out,
                  delta_r,
                  bounds,
                  nindent=0):
    """
    Recursively search for objects in the JSON file 
    containing a 'ingredients' field.  
    The 'ingredients' should map to a dictionary of molecule_tyle_names.
    For each of them, define a moltemplate molecule type, followed by
    a list of \"new\" commands which instantiate a copy of that molecule
    at different positions and orientations.
    """

    if not isinstance(tree, dict):
        return
    if 'ingredients' in tree:
        ConvertMolecules(tree['ingredients'],
                         file_out,
                         delta_r,
                         bounds,
                         nindent)
    else:
        for object_name in tree:
            if object_name == 'recipe':
                continue
            file_out.write(nindent*'  '+object_name + ' {\n')
            ConvertSystem(tree[object_name],
                          file_out,
                          delta_r,
                          bounds,
                          nindent+1)
            file_out.write(nindent*'  '+'}  # endo of \"'+object_name+'\" definition\n\n')
            file_out.write('\n' + 
                           nindent*'  '+object_name + '_instance = new ' + object_name + '\n' +
                           '\n' +
                           '\n')

            #when debugging, uncomment the next line:
            #break





def ConvertCellPACK(file_in,        # typically sys.stdin
                    filename_in,    # optional file name
                    file_out,       # typically sys.stdout
                    filename_out,   # optional file name
                    out_obj_name,   # optional moltemplate object name
                    delta_r,        # radial resolution
                    pairstyle,      # which LAMMPS pair_style do you want?
                    pairstyle2docs, # documentation for these pair styles
                    pairstyle2args, # required arguments for these pair styles
                    epsilon,        # Lennard-Jones parameter
                    debye):         # Debyle length (if applicable)
    """
    Read a JSON file created by CellPACK and 
    convert it to MOLTEMPLATE (LT) format.
    """
    tree = json.load(file_in)
    # depth-first-search the JSON tree looking
    # for any nodes containing the string 'ingredients'
    # convert them to molecule definitions and
    # lists of instances (copies)
    if out_obj_name != '':
        file_out.write('# This file defines a type of molecular object (\"'+out_obj_name+'\") contaning\n'
                       '# the entire system from the JSON output of CellPACK.  Later on, before you\n'
                       '# use moltemplate.sh, you must create a file (eg \"system.lt\") containing:\n'
                       '# import \"'+filename_out+'\"\n'
                       '# system = new '+out_obj_name+'\n'
                       '#  OR, if you want multiple copies of your system, replace the line above with:\n'
                       '# copy_1 = new '+out_obj_name+'.rot(theta1,ax1,ay1,az1).move(x1,y1,z1)\n'
                       '# copy_2 = new '+out_obj_name+'.rot(theta2,ax2,ay2,az2).move(x2,y2,z2)\n'
                       '# copy_3 = new '+out_obj_name+'.rot(theta3,ax3,ay3,az3).move(x3,y3,z3)\n'
                       '#     :\n'
                       '# Then run moltemplate.sh on the \"system.lt\" file using:\n'
                       '#\n'
                       '# moltemplate.sh -nocheck '+filename_out+'\n')
                       #'#      (\"moltemplate.sh system.lt\" works too but takes longer)\n')
    else:
        file_out.write('# Later, you can run moltemplate on this file using\n'
                       '# moltemplate.sh -nocheck '+filename_out+'\n')
                       #'#       (moltemplate.sh '+filename_out+' works too but takes longer)\n')
    #file_out.write('# On a large system, moltemplate.sh can up to 2 hours to complete.\n'





    file_out.write('\n'
                   '\n')
    nindent = 0
    if out_obj_name != '':
        file_out.write(out_obj_name + ' {\n')
        nindent = 1


    file_out.write('\n\n'
                   'ForceField {\n'
                   '\n'
                   '   # The \"ForceField\" object defines a list of atom types\n'
                   '   # and the forces between them.  These atom types and forces\n'
                   '   # are shared by all of the molecule types in this simulation.\n'
                   '\n')




    # Figure out what kinds of particles we need
    # These will be referenced below when we
    # create a "ForceField" object which describes them.

    ir_needed = set([]) # Which particle radii are needed in this simulation?
                        # LAMMPS does not allow every particle in the sim to
                        # have a unique radius.  The number of particle types
                        # cannot exceed roughly 10^3 (for efficient execution)
                        # Hence, the radii of each particle is quantized by
                        # rounding it to the nearest multiple of "delta_r".
                        # All particles whose radii fall into this bin
                        # are represented by the same particle type.
                        # The "ir_needed" keeps track of which bins have been
                        # visited.  Later, we will loop through all the bins
                        # in "ir_needed", and define particle types for the
                        # corresponding particles that fall into that bin
                        # (that range of radii sizes), and choose appropriate
                        # force-field parameters (Lennard-Jones parameters)
                        # for those particle types.

    RadiiNeeded(tree, ir_needed, delta_r)
    assert(len(ir_needed) > 0)



    rmax = max(ir_needed) * delta_r
    rcut = rmax * (2.0**(1.0/6))
    bounds = [[0.0,0.0,0.0],    #Box big enough to enclose all the particles
              [-1.0,-1.0,-1.0]] #[[xmin,ymin,zmin],[xmax,ymax,zmax]]
    pairstyle2args['lj/cut'] = str(rcut)
    pairstyle2args['lj/cut/coul/debye'] = str(debye)+' '+str(rcut)+' '+str(rcut)
    pairstyle2args['lj/cut/coul/cut'] = str(rcut)
    pairstyle2args['lj/cut/coul/long'] = str(rcut)
    pairstyle2args['lj/class2/coul/long'] = str(rcut)
    pairstyle2args['lj/class2/coul/cut'] = str(rcut)
    #pair_mixing_style = 'sixthpower tail yes'
    pair_mixing_style = 'arithmetic'
    special_bonds_command = 'special_bonds lj/coul 0.0 0.0 1.0'



    file_out.write('   write_once("In Settings") {\n')
    for iradius in sorted(ir_needed):
        rcut = 2 * iradius * delta_r   #(don't forget the 2)
        r = rcut / (2.0**(1.0/6))
        file_out.write('     pair_coeff ' +
                       '@atom:A' + str(iradius) + ' ' +
                       '@atom:A' + str(iradius) + ' ' +
                       pairstyle + ' ' +
                       str(epsilon) + ' ' +
                       str(r) + ' ' +
                       str(rcut) + '\n')
    file_out.write('  }  #end of "In Settings" (pair_coeff commands)\n'
                   '\n')
    default_mass = 1.0   # Users can override this later
    file_out.write('  # Last I checked, LAMMPS still requires that every atom has its\n'
                   '  # mass defined in the DATA file, even if it is irrelevant.\n'
                   '  # Take care of that detail below.\n')
    file_out.write('\n'
                   '   write_once("Data Masses") {\n')

    for iradius in sorted(ir_needed):
        rcut = iradius * delta_r
        r = rcut / (2.0**(1.0/6))
        file_out.write('    ' +
                       '@atom:A' + str(iradius) +' '+ str(default_mass) +'\n')
    file_out.write('  }  # end of "Data Masses"\n')
    file_out.write('\n\n'
                   '  # At some point we must specify what -kind- of force fields we want to use\n'
                   '  # We must also specify how we want to measure distances and energies (units)\n'
                   '  # as well as how to represent the particles in the simulation (atom_style)\n'
                   '\n'
                   '  write_once("In Init") {\n'
                   '\n'
                   '    atom_style full #(default atom style)\n'
                   '\n'
                   '    units lj        #(this means the units can be customized by the user (us))\n'
                   '\n'
                   '    pair_style hybrid ' + pairstyle + ' ' +
                   pairstyle2args[pairstyle] + '\n')
    file_out.write('    # For details, see\n' +
                   '    # ' + pairstyle2docs[pairstyle] + '\n')

    file_out.write('    # Use ordinary Lorenz-Berthelot mixing rules.\n'
                   '    # (This means the effective minimal distance\n'
                   '    #  between unlike particles is the sum of their radii,\n'
                   '    #  as it would be if they were hard spheres)\n'
                   '\n'
                   '    pair_modify mix ' + pair_mixing_style + '\n')
    file_out.write('\n'
                   '    # The next line is probably not relevant, since we\n'
                   '    # do not have any bonds in our system.\n'
                   '    ' + special_bonds_command + '\n')
    file_out.write('\n'
                   '    # If I am not mistaken, the next two lines help speed up the computation\n'
                   '    # of neighbor lists.  They are necessary because of the wide diversity\n'
                   '    # of particle sizes in the simulation.  (Most molecular dynamics software\n'
                   '    # was optimized for particles of similar size.)\n'
                   '\n'
                   '    comm_modify mode multi  # Needed since we have many different particle sizes\n'
                   '    neighbor 3.0 multi      # Adjust this number later to improve efficiency\n'
                   '\n')
    file_out.write('  }  # finished selecting force field styles\n')

    file_out.write('\n'
                   '\n'
                   '  # Optional: Create a file to help make display in VMD prettier\n'
                   '\n'
                   '  write_once("vmd_commands.tcl") {\n')

    for iradius in sorted(ir_needed):
        r = iradius * delta_r
        atype_name = '@{atom:A' + str(iradius) + '}'
        file_out.write('    set sel [atomselect top "type '+atype_name+'"]\n'
                       '    \$sel set radius '+str(r)+'\n')

    file_out.write('  }  # end of "vmd_commands.tcl"\n\n')

    file_out.write('\n'
                   '\n'
                   '  ### Note: Use a rigid-body integrator to keep the particles\n'
                   '  ###      in each molecule from moving relative to eachother:\n'
                   #'\n'
                   #'  #write_once("In Settings") {\n'
                   '  # group gRigid @atom:.../ForceField/*\n'
                   '  # http://lammps.sandia.gov/doc/group.html\n'
                   '  #\n'
                   '  # fix fxRigid gRigid rigid molecule\n'
                   '  # http://lammps.sandia.gov/doc/fix_rigid.html\n'
                   '  #\n'
                   '  ### This next line greatly increases the speed of the simulation:\n'
                   '  ### No need to calculate forces between particles in the same rigid molecule\n'
                   '  # neigh_modify exclude molecule/intra gRigid\n'
                   '  # http://lammps.sandia.gov/doc/neigh_modify.html\n'
                   #'  #}\n'
                   '\n')
    file_out.write('}  # end of the "ForceField" object definition\n'
                   '\n\n\n\n\n')


    ConvertSystem(tree['compartments'],
                  file_out,
                  delta_r,
                  bounds,
                  nindent)

    if out_obj_name != '':
        file_out.write('}  # end of \"'+ out_obj_name + '\" object definition\n')

    

    # Print the simulation boundary conditions

    for d in range(0, 3):      # Simulation particles are not point-like
        bounds[0][d] -= rcut   # add extra space to the boundary box to 
        bounds[1][d] += rcut   # accomodate particles of finite size (rcut)

    file_out.write('\n\n'
                   '# Simulation boundaries:\n'
                   '\n'
                   'write_once("Data Boundary") {\n'
                   '  '+str(bounds[0][0])+' '+str(bounds[1][0])+' xlo xhi\n'
                   '  '+str(bounds[0][1])+' '+str(bounds[1][1])+' ylo yhi\n'
                   '  '+str(bounds[0][2])+' '+str(bounds[1][2])+' zlo zhi\n'
                   '}\n\n')








def main():
    try:
        sys.stderr.write(g_program_name + ", version " +
                         __version__ + ", " + __date__ + "\n")

        if sys.version < '2.7':
            from ordereddict import OrderedDict
        else:
            from collections import OrderedDict 
    
        if sys.version > '3':
            import io
        else:
            import cStringIO
    
        # defaults:
        out_obj_name = ''
        type_subset = set([])
        filename_in = ''
        filename_out = 'THIS_FILE'
        file_in = sys.stdin
        file_out = sys.stdout

        delta_r = 0.1     # resolution of particle radii

        # --- Units ---
        # The following parameters depend on what units you are using.
        # By default, this program assumes distances are in nm,
        # energies are in kCal/Mole,
        # Note: When running LAMMPS, temperatures in the fix_langevin or
        #       fix_nvt commands must be specified by the user in units of
        #       energy (ie kCal/mole), not K.   (In these units,
        #       k_B*temperature = 0.5961621, assuming temperature = 300K)
        debye = 1.0       # the Debye length(only relevent for some pair_styles)
        kB = 0.001987207    # Default Boltzmann's constant ((kCal/Mole)/degreeK)
        temperature = 300.0 # Default temperature (in K)
        epsilon = kB*temperature  # The Lennard-Jones "epsilon" parameter has
                                  # units of energy and should be approximately
                                  # equal to the value of k_B*temperature in 
                                  # whatever units you are using.
        pairstyle = 'lj/cut'
        pairstyle2docs = {}
        pairstyle2args = defaultdict(str)
        pairstyle2docs['lj/cut'] = 'http://lammps.sandia.gov/doc/pair_lj.html'
        pairstyle2docs['lj/cut/coul/debye'] = 'http://lammps.sandia.gov/doc/pair_lj.html'
        pairstyle2docs['lj/cut/coul/cut'] = 'http://lammps.sandia.gov/doc/pair_lj.html'
        pairstyle2docs['lj/cut/coul/long'] = 'http://lammps.sandia.gov/doc/pair_lj.html'
        pairstyle2docs['lj/class2/coul/long'] = 'http://lammps.sandia.gov/doc/pair_class2.html'
        pairstyle2docs['lj/class2/coul/cut'] = 'http://lammps.sandia.gov/doc/pair_class2.html'

        argv = [arg for arg in sys.argv]
    
        i = 1
    
        while i < len(argv):
    
            if argv[i] == '-name':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by the name of the\n'
                                     '       moltemplate object you wish to create.\n')
                out_obj_name = argv[i + 1]
                del argv[i:i + 2]

            elif argv[i] == '-pairstyle':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by a LAMMPS pair_style.\n')
                pairstyle = argv[i + 1]
                if not pairstyle in pairstyle2docs:
                    raise InputError('Error: Invalid '+argv[i]+' argument.  Available choices are:\n'
                                     '\n'.join([ps for ps in pairstyle2docs])+'\n')
                del argv[i:i + 2]

            elif argv[i] in ('-deltaR', '-deltar', '-delta-r'):
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by a number.\n')
                delta_r = float(argv[i + 1])
                del argv[i:i + 2]

            elif argv[i] == '-epsilon':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by a number.\n')
                epsilon = float(argv[i + 1])
                del argv[i:i + 2]

            elif argv[i] == '-debye':
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by a number.\n')
                debye = float(argv[i + 1])
                del argv[i:i + 2]

            elif argv[i] in ('-in', '-file'):
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by\n'
                                     '       the name of a JSON file created by CellPACK.\n')
                filename_in = argv[i + 1]
                try:
                    file_in = open(filename_in, 'r')
                except IOError:
                    sys.stderr.write('Error: Unable to open file\n'
                                     '       \"' + filename_in + '\"\n'
                                     '       for reading.\n')
                    sys.exit(1)
                del argv[i:i + 2]

            elif argv[i] in ('-out'):
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by the name of a moltemplate file\n')
                filename_out = argv[i + 1]
                try:
                    file_out = open(filename_out, 'w')
                except IOError:
                    sys.stderr.write('Error: Unable to open file\n'
                                     '       \"' + filename_out + '\"\n'
                                     '       for writing.\n')
                    sys.exit(1)
                del argv[i:i + 2]

            elif argv[i] in ('-url', '-in-url'):
                import urllib2
                if i + 1 >= len(argv):
                    raise InputError('Error: ' + argv[i] + ' flag should be followed by a URL\n'
                                    '       pointing to a JSON file created by CellPACK.\n')

                url = argv[i + 1]
                try:
                    request = urllib2.Request(url)
                    file_in = urllib2.urlopen(request)
                except urllib2.URLError:
                    sys.stdout.write("Error: Unable to open link:\n" + url + "\n")
                    sys.exit(1)
                del argv[i:i + 2]
    
            elif argv[i] in ('-help', '--help', '-?', '--?'):
                sys.stderr.write(doc_msg)
                sys.exit(0)
                del argv[i:i + 1]
    
            else:
                i += 1
    
        if len(argv) != 1:
            raise InputError('Error: Unrecongized arguments: ' + ' '.join(argv[1:]) +
                             '\n\n' + doc_msg)



        file_out.write('# This file was generated automatically using:\n'+
                       '#    '+g_program_name+' '+' '.join(sys.argv[1:])+'\n')


        ConvertCellPACK(file_in,
                        filename_in,
                        file_out,
                        filename_out,
                        out_obj_name,
                        delta_r,
                        pairstyle,
                        pairstyle2docs,
                        pairstyle2args,
                        epsilon,
                        debye)


    except InputError as err:
        sys.stderr.write('\n\n' + str(err) + '\n')
        sys.exit(1)

if __name__ == '__main__':
    main()
