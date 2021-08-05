"""
Utility functions implemented here are clearly defined unit operations. They may only be called once in the code, but
are defined separately for cleanliness and legibility.
"""

import os
import re
import sys
import copy
import math
import numpy
#import pytraj
#import mdtraj
#import parmed
import argparse
#from paprika.build.system import TLeap as tleap
import fileinput
#import dill as pickle   # I think this is kosher!
#from simtk.openmm.app import *
#from simtk.openmm import *
#from simtk.unit import *
#from isee.initialize_charges import set_charges

# Two different ways to import tleap depending on I think paprika version
# todo: figure out why this is necessary and how to make ModuleNotFoundError accessible
try:
    from paprika.build.system import TLeap as tleap
except ModuleNotFoundError:
    from paprika import tleap



def build_peptide(thread, settings):
    """
    Build peptide according to amino acid sequence in Settings.sequence using Amber's TLEaP

    Parameters
    ----------
    thread : Thread
        Thread on which to operate
    settings : argparse.Namespace
        Settings namespace object


    Returns
    -------
    tleap init_coords or None?

    """
    # todo: Call to TLEaP. Determine where/how files will be outputted. Either return file name of .pdb file generated
    # todo: or directly add to settings.init_coord(?) w/in this function
    # write sequence in tleap format, e.g. ['ACE', 'LYS', 'NME'] --> '{ ACE LYS NME }
    sequence = '{ ' + " ".join(thread.peptides[0])) + ' }' # todo: how to index into this?

    try:
        system = tleap()
    except TypeError:
        system = tleap.System()
    system.pbc_type = None
    system.neutralize = False
    system.output_path = settings.working_directory # todo: check on this
    system.output_prefix = thread.name + '_tleap' # todo: will need to specify which peptide also
    system.template_lines = [
        'source oldff/leaprc.ff99SBildn', # todo: update ff, not completely necessary
        thread.name ' = sequence' + sequence, # todo: sequence = thread.peptides[i]
        'savepdb ' thread.name str(thread.peptide[0]) + '.pdb', # todo: decide on file naming
        'quit',
    ]
    return None

def edit_pdb(thread, settings):
    """
    Edits .pdb file from TLEaP to make compatible with Gromacs method pdb2gmx
    Parameters
    ----------
    thread : Thread
        Thread on which to operate
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------

    """
    pdb = thread.name + '_tleap' # todo: check naming based on tleap output
    # todo: test this method further
    previous_line = ''
    h_list = ['HA', 'HB', 'HD', 'HE', 'HG']
    h_prev = 0
    residues = {}

    with open(pdb_in, 'r') as pdb_in, open(pdb_out, 'w') as pdb_out:
        # head = [next(pdb_in) for x in range(len(pdb_in))[:-2]]
        for line in pdb_in:
            line_split = re.split(r'(\s+)', line)

            # skip last two lines
            if line_split[0] == 'TER' or line_split[0] == 'END':
                continue

            atom = line_split[4]
            AA = line_split[6]
            res = line_split[8]

            # change histidine naming
            if AA == 'HIE' or AA == 'HIP' or AA == 'HID': # todo: determine how best to deal w/ protonation
                line_split[6] = 'HIS'

            # special case for  SER, CYS
            if (AA == 'CYS' or AA == 'SER') and atom == 'HG':
                atom = 'HG1'
                # reformat white space
                line_split[5] = ' ' * (len(line_split[5]) - 1)
            else:
                # if hydrogren number does not begin with 1
                if len(atom) >= 3 and atom[0:2] in h_list:
                    h_num = int(atom[-1])
                    if h_num != 1:
                        h_mod = str(h_prev + 1)
                        if atom[:-1] + h_mod not in residues[res]:
                            atom = atom[:-1] + h_mod
                        h_prev = int(h_mod)
                    else:
                        h_prev = int(h_num)
                else:  # reset h_prev
                    h_prev = 0
            # reassign line values
            line_split[4] = atom

            # add to dictionary to track used atom names
            residues.setdefault(res, []).append(atom)

            pdb_out.write(''.join(line_split))


def pdb2gmx(thread, settings):
    """

    Parameters
    ----------
    thread : Thread
        Thread on which to operate
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------

    """
    # todo: check that subprocess is the best way of handling; determine location of inputs
    commandline_arg = 'gmx_mpi pdb2gmx -f {} -o {} -p {} -ter < {}'.format(pdb_file, gro_file, protein_ff, user_inp)
    subprocess.run(commandline_arg, shell=True)
    # todo: need to write and store user input file in inputs

def clean_peptide_ff(thread, settings):
    """
    Moves Gromacs generated peptide topology from topol.top to {peptide}.itp. Writes new topol.top to contain
    information on the force fields to include and eventually to add

    Parameters
    ----------
    thread
    settings

    Returns
    -------

    """

def add_to_topology(thread, settings):
    """
    Writes to base topol.top to reference force fields to be included (user selection,

    Parameters
    ----------
    thread
    settings

    Returns
    -------

    """

def grow_box(thread, settings):
    # todo: testing needed
    # read file and set new box size to dimenstions (x y z) nm
    lines = open(gromacs_file, 'r').readlines()
    new_box = ('   {}   {}   {}'.format(x, y, z))
    lines[-1] = new_box

    # reopen and write to file
    open(gromacs_file, 'w').writelines(lines)

def center_peptide(thread, settings):
    # todo: testing needed
    commandline_arg = 'gmx_mpi editconf -f {} -o {} -center {} {} {}'.format(gro_file, output_file, x, y, z)
    subprocess.run(commandline_arg, shell=True)

def solvate(thread, settings):
    # todo: testing needed
    commandline_arg = 'gmx_mpi solvate -cp {} -cs spc216.gro -o {} -p {}'.format(gro_file, output_file, topol_file)
    subprocess.run(commandline_arg, shell=True)

def grompp_ion_runfile(thread, settings):
    # todo testing needed
    commandline_arg = 'gmx_mpi grompp -f ion.mdp -c {} -p topol.top -o {} -maxwarn 4'.format(gro_file, tpr_file)
    subprocess.run(commandline_arg, shell=True)

def genion(thread, settings):
    # todo testing needed
    commandline_arg = 'gmx_mpi genion -s {} -o {} -p topol.top -pname NA -nname CL -neutral -np 5 < {}'.format(tpr_file, output_file, user_inp)
    subprocess.run(commandline_arg, shell=True)
    # todo: need user input file

def translate(thread, settings):
    # todo: testing needed
    commandline_arg = 'gmx_mpi editconf -f {} -o {} -translate {} {} {}'.format(gro_file, output_file, x, y, z)
    subprocess.run(commandline_arg, shell=True)

def convert_coord(thread, settings):
    # todo: testing needed
    # converting .pdb to .gro or .gro to .pdb
    commandline_arg = 'gmx_mpi editconf -f {} -o {}'.format(file_in, file_out)
    subprocess.run(commandline_arg, shell=True)

def create_index(thread, settings):
    # todo: testing needed. Input file 'q' needed
    commandline = 'gmx_mpi make_ndx -f {} -o {}'.format(gro_file, index_file)
    subprocess.run(commandline_arg, shell=True)

def combine_index(thread, settings):
    # todo: more pythonic way of doing this than cat!
    pass



