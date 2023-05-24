#!/usr/bin/env python3

import sys
import numpy as np

"""
    Open and read log files of frequency calculations.
    Current options: g09 (GAUSSIAN09)
"""

def catch_hp_g09(LogFile):
    """
        Determine if the g09 log file corresponds to a High 
        Precision calculation (bool)
    """

    with open(LogFile) as f:

        for line in f.readlines():
            if "HPMODES" in line or "hpmodes" in line:
                return True
        return False

def catch_N_atoms_g09(LogFile):

    """ 
        Catch number of atoms of the molecule
    """
    
    with open(LogFile) as f:

        for line_number, line in enumerate(f):
            if "Input orientation:" in line:
                geometryLine0 = line_number+5
            
            if "Distance matrix" in line:
                geometryFinalLine = line_number-1
                break

    N_atoms = geometryFinalLine - geometryLine0
    return N_atoms

def catch_rot_symmetry_number(LogFile):
    """
        Catch rotational symmetry number
    """

    with open(LogFile) as f:

        for line_number, line in enumerate(f):
            if "Rotational symmetry number" in line:
                rot_line = line

                break

        rot_number = int(str(rot_line.split()[-1])[:-1])
        return rot_number

def catch_atom_list_g09(LogFile,N_atoms):
 
    """ 
        Catch atom list (list)
    """

    atom_list = []

    with open(LogFile) as f:

        for line_number, line in enumerate(f):
            if "Charge =" in line:
                atom_list_line0 = line_number+1
            
    with open(LogFile) as f:

        for line_number, line in enumerate(f):

            if line_number in range(atom_list_line0,atom_list_line0+N_atoms):
                atom_list.append(line.split()[0])    

    return atom_list

def catch_geom_g09(LogFile):

    """ 
        Catch Geometry in angs (ndarray)
    """

    N_atoms = catch_N_atoms_g09(LogFile)
    geom = np.ndarray(shape=(N_atoms,3),dtype=np.float32)
    
    with open(LogFile) as f:

        for line_number, line in enumerate(f):
            if "Input orientation:" in line:
                geometryLine0 = line_number+5

            if "Distance matrix" in line:
                geometryFinalLine = line_number-1

    with open(LogFile) as f:

        for line_number, line in enumerate(f):
            if line_number in range(geometryLine0,geometryFinalLine):

                geom[line_number-geometryLine0] = [float(line.split()[3]),
                                                   float(line.split()[4]),
                                                   float(line.split()[5])]

    return geom

def catch_freq_g09(LogFile,hp=True):

    """
        Catch frequencies (list) (cm -1)
        The keyword hp search for frequencies in High Precision calculations (default)
        ---> g09 input: freq=(hpmodes)
    """

    freq = []

    if hp:
        with open(LogFile) as f:

            for line in f.readlines():
                if "Frequencies ---" in line:

                    for i in range(2,len(line.split())):
                        freq.append(float(line.split()[i]))        

    else:
        with open(LogFile) as f:

            for line in f.readlines():
                if "Frequencies -- " in line:

                    for i in range(2,len(line.split())):
                        freq.append(float(line.split()[i]))  

    return freq

def catch_rot_temp(LogFile):
    """
        Catch rotational temperatures
    """

    with open(LogFile) as f:
        
        for line_number, line in enumerate(f):
            if "Rotational temperatures (Kelvin)" in line:
                Ix = float(line.split()[-3])
                Iy = float(line.split()[-2])
                Iz = float(line.split()[-1])
                break

            if "Rotational temperature (Kelvin)" in line:
                Ix = float(line.split()[-1])
                Iy = float(line.split()[-1])
                Iz = float(line.split()[-1])
                break
                
        return Ix, Iy, Iz

def catch_spin_mult(LogFile):
    """
        Catch spin multiplicity
    """

    with open(LogFile) as f:

        for line_number, line in enumerate(f):
            if "Multiplicity =" in line:
                spin_line = line

                break

        spin_mult = int(spin_line.split()[-1])
        return spin_mult

def catch_mol_mass(LogFile):
    """
        Catch total molecular mass
    """

    with open(LogFile) as f:

        for line_number, line in enumerate(f):
            if "Molecular mass:" in line:
                mass_line = line

                break

        mass = float(mass_line.split()[-2])
        return mass

def catch_red_mass_g09(LogFile,hp=True):

    """
        Catch reduced masses (list) (a.m.u)
        The keyword hp search for frequencies in High Precision calculations (default)
        ---> g09 input: freq=(hpmodes)
    """

    red_mass = []

    if hp:
        with open(LogFile) as f:

            for line in f.readlines():
                if "Reduced masses ---" in line:

                    for i in range(3,len(line.split())):
                        red_mass.append(float(line.split()[i]))        

    else:
        with open(LogFile) as f:

            for line in f.readlines():
                if "Red. masses -- " in line:

                    for i in range(3,len(line.split())):
                        red_mass.append(float(line.split()[i]))  

    return red_mass

def catch_ir_intensity_g09(LogFile,hp=True):

    """
        Catch Infra-Red intensities (list) (km/mol)
        The keyword hp search for frequencies in High Precision calculations (default)
        ---> g09 input: freq=(hpmodes)
    """

    ir_intensity = []

    if hp:
        with open(LogFile) as f:

            for line in f.readlines():
                if "IR Intensities ---" in line:

                    for i in range(3,len(line.split())):
                        ir_intensity.append(float(line.split()[i]))        

    else:
        with open(LogFile) as f:

            for line in f.readlines():
                if "IR Inten    -- " in line:

                    for i in range(3,len(line.split())):
                        ir_intensity.append(float(line.split()[i]))  

    return ir_intensity

def catch_displacement_vectors_g09(LogFile,hp=True):

    """
        Catch sequence of displacement vectors (ndarray ---> all normal coordinates)
        The format of the list is 3N-6 (or 3N-5 for linear molecules) ndarrays of 
        N_atoms x 3 dimensions each.
        The keyword hp search for frequencies in High Precision calculations (default)
        ---> g09 input: freq=(hpmodes)
    """

    N_atoms = catch_N_atoms_g09(LogFile)
    N_modes = len(catch_freq_g09(LogFile,hp))
    displacement_vectors = np.ndarray(shape=(N_modes,N_atoms,3),dtype=np.float32)

    if hp:

        with open(LogFile) as f:

            for line_number, line in enumerate(f):
                if "Coord Atom Element:" in line:
                    mode_labels_line = line_number-6
                    displacement_Line0 = line_number+1

                    # New loop #
                    with open(LogFile) as newf:

                        for newline_number, newline in enumerate(newf):
                            if newline_number == mode_labels_line:
                                mode_labels = newline.split()
                            
                            if newline_number in range(displacement_Line0,displacement_Line0+3*N_atoms):
                                for k in mode_labels:
                                    displacement_vectors[int(k)-1,
                                    int(np.floor((newline_number-displacement_Line0)/3)),
                                    int((newline_number-displacement_Line0)%3)] = float(newline.split()[3+int(k)-int(mode_labels[0])])
                    # End of new loop #
        
    else:
        
        with open(LogFile) as f:

            for line_number, line in enumerate(f):
                if "Atom  AN" in line:
                    mode_labels_line = line_number-6
                    displacement_Line0 = line_number+1

                    # New loop #
                    with open(LogFile) as newf:

                        for newline_number, newline in enumerate(newf):
                            if newline_number == mode_labels_line:
                                mode_labels = newline.split()                            

                            if newline_number in range(displacement_Line0,displacement_Line0+N_atoms):
                                for k in mode_labels:
                                    for j in range(3):
                                        displacement_vectors[int(k)-1,
                                        int(newline_number-displacement_Line0),
                                        j] = float(newline.split()[2+3*(int(k)-int(mode_labels[0]))+j])
                    # End of new loop #

    return displacement_vectors
