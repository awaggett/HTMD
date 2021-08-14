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
import subprocess
#from paprika.build.system import TLeap as tleap
import fileinput
#import dill as pickle   # I think this is kosher!
#from simtk.openmm.app import *
#from simtk.openmm import *
#from simtk.unit import *
#from htmd_.initialize_charges import set_charges

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

    # or for bash script purposes: input seq = thread.peptides[thread.peptide]
    # sequence = '{ ' + " ".join(seq) + ' }'

    sequence = '{ ' + " ".join(thread.peptides[thread.peptide]) + ' }'

    try:
        system = tleap()
    except TypeError:
        system = tleap.System()
    short_name = thread.name + '_' + thread.peptide
    long_name = short_name + '_tleap.pdb'
    system.pbc_type = None
    system.neutralize = False
    system.output_path = settings.working_directory # todo: check on this
    system.output_prefix = thread.name + '_' + thread.peptide + '_tleap' # todo: will need to specify which peptide also
    system.template_lines = [
        'source oldff/leaprc.ff99SBildn', # todo: update ff, not completely necessary
        short_name + ' = sequence ' + sequence, # todo: sequence = thread.peptides[i]
        'savepdb ' + short_name + ' ' + long_name, # todo: decide on file naming
        'quit',
    ]

    # add output coordinate file to thread.history.coords
    thread.history.coords.append(long_name)
    return long_name

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
    # Extract pdb from thread.history.coords
    pdb_in = thread.history.coords[-1]
    # todo: test this method further
    pdb_out = pdb_in.split('.')[0] + '_mod.pdb'
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

        # todo: add coord file to thread.history.coords
        thread.history.coords.append(pdb_out)
    return pdb_out

def get_input(input_file):
    if os.path.exists(settings.path_to_inout_files + '/' + input_file):
        return input_file # todo: do I need the full path here?

def get_template(templ):
    if os.path.exists(settings.path_to_templates + '/' + templ):
        return templ # todo: do I need the full path here?

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
    pdb_file = thread.history.coords[-1]
    gro_file = thread.name + '_' + thread.peptide + '_init.gro'
    protein_ff = settings.force_field.split('.')[0]
    topol_file = thread.name + '_' + thread.peptide + '.top'

    commandline_arg = 'echo 3 4 | gmx_mpi pdb2gmx -f {} -o {} -ff {} -p {} -water none -ter '.format(pdb_file, gro_file, protein_ff, topol_file)
    subprocess.run(commandline_arg, shell=True)

    thread.history.coords.append(gro_file)
    thread.history.tops.append(topol_file)
    # todo: can maybe also specify naming of topology file here to add to thread history (and posre to remove?)
    return gro_file

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
    # todo: remove head and tail from generated topol.top and move to .itp
    # todo: will I need to extract the #include statements in case naming of file paths is different between force fields?
    # create new itp file for protein:q

    protein_ff_in = thread.history.tops[-1]
    protein_ff_out = thread.name + '_' + thread.peptide + '.itp'
    ff_out = open(protein_ff_out, 'w')

    # remove unnecesary lines from gromacs generated topology file and write to .itp
    with open(protein_ff_in, 'r') as ff_in:

        start_writing = False
        end_writing = False
        for line in ff_in:
            if '[ moleculetype ]' in line:
                start_writing = True
                ff_out.write(line)
            elif 'Position restraint' in line: # todo: naming is different I think - need to check this
                end_writing = True
            elif start_writing == True and end_writing == False:
                ff_out.write(line)

    # remove gromacs generated topology file and posre file # todo: is there a way to avoid writing these?
    # os.remove('posre.itp')
    # os.remove(protein_ff_in)


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
    # need template topology file here - this is necessary because only way in gromacs to later combine topologies
    # create new topology file
    top_out = open(topol_file, 'w')
    topol_base = get_template('topol.top')

    # add protein ff to topology file
    with open(topol_base, 'r') as top_in:

        for line in top_in:
            line_split = re.split(r'(\s+)', line)
            if line_split[0] == '#protein_ff':
                top_out.write('#include "./{}.itp"'.format(base_name))
            elif line_split[0] == 'sys_name':
                top_out.write('{} on Quartz slab'.format(base_name))
            elif line_split[0] == 'protein_name':
                top_out.write('{}       1'.format(base_name))
            else:
                top_out.write(line)

#def grow_box(thread, settings):
#    # todo: this can be accomplished with centering
#    # read file and set new box size to dimenstions (x y z) nm
#    lines = open(gromacs_file, 'r').readlines()
#    new_box = ('   {}   {}   {}'.format(settings.peptide_box_dim, settings.peptide_box_dim, settings.peptide_box_dim))
#    lines[-1] = new_box
#    # todo: check if new box is meant to be a tuple...
#
#    # reopen and write to file
#    open(gromacs_file, 'w').writelines(lines)

def get_surface_size(settings):
    surface = settings.surface_coord
    dims = open(surface, 'r').readlines()[-1].split()
    return dims[0], dims[1]
    # todo: check this returns the correct thing

def center_peptide(thread, settings):
    # todo: testing needed !!! -center can also handle growing box
    gro_file = thread.history.coords[-1]
    output_file = thread.name + '_' + thread.peptide + '_' + thread.system + '_center.gro' # todo: decide on naming convention
    dim = settings.peptide_box_dim
    commandline_arg = 'gmx_mpi editconf -f {} -o {} -box {} '.format(gro_file, output_file, dim)
    subprocess.run(commandline_arg, shell=True)
    thread.history.coords.append(output_file)

def solvate(thread, settings):
    # todo: testing needed
    gro_file = thread.history.coords[-1]
    output_file = thread.name + '_' + thread.peptide + '_' + thread.system + '_solvate.gro' # todo: decide on naming convention
    topol_file = thread.history.tops[-1]
    commandline_arg = 'gmx_mpi solvate -cp {} -cs spc216.gro -o {} -p {}'.format(gro_file, output_file, topol_file)
    subprocess.run(commandline_arg, shell=True)
    # todo: should keep same naming of topology. will just need unique naming for thread.name + thread.peptide + thread.system
    thread.history.coords.append(output_file)

def grompp_ion_runfile(thread, settings):
    # todo testing needed
    gro_file = thread.history.coords[-1]
    tpr_file = thread.name + '_' + thread.name + '_' + thread.system + '_ion.tpr' # todo: check naming!
    commandline_arg = 'gmx_mpi grompp -f ion.mdp -c {} -p topol.top -o {} -maxwarn 4'.format(gro_file, tpr_file)
    subprocess.run(commandline_arg, shell=True)
    thread.history.runfiles.append(tpr_file)

def get_system_charge(thread, settings):
    # maybe there is a gromacs command for this?
    pass

def genion(thread, settings):
    # todo testing needed
    tpr_file = thread.history.runfiles[-1]
    output_file = thread.name + '_' + thread.name + '_' + thread.system + '_ion.gro' # todo: check!
    input_file = get_input('genion_input.txt')
    charge = get_system_charge(thread, settings)
    if charge > 0:
        charge_mod = '-np 5'
    else:
        charge_mod = '-nn 5'
    commandline_arg = 'gmx_mpi genion -s {} -o {} -p topol.top -pname NA -nname CL -neutral {} < {}'.format(tpr_file, output_file, charge_mod, input_file)
    subprocess.run(commandline_arg, shell=True)
    thread.history.coords.append(output_file)
    # todo need to decide on number np or nn added based on system charge

# todo: actually, I am going to keep this in batch job
def grompp_run(thread, settings):
    gro_file = thread.history.coords[-1]
    top_file = thread.history.tops[-1]
    output_file = thread.name + thread.suffix # todo: check!

    if thread.current_type == 'peptide_em':
        input_file =  get_input('em_pep.mdp')
        commandline = 'gmx_mpi grompp -f {} -c {} -p {} -o {} -maxwarn 3'.format(input_file, gro_file, top_file, output_file)
    elif thread.current_type == 'peptide_npt':
        input_file =  get_input('npt_pep.mdp')
        commandline = 'gmx_mpi grompp -f {} -c {} -p {} -o {} -maxwarn 3'.format(input_file, gro_file, top_file, output_file)
    elif thread.current_type == 'peptide_nvt':
        input_file =  get_input('nvt_pep.mdp')
        commandline = 'gmx_mpi grompp -f {} -c {} -p {} -o {} -maxwarn 3'.format(input_file, gro_file, top_file, output_file)
    elif thread.current_type == 'system_em':
        input_file = get_input('em_sys.mdp')
        index_file = thread.history.indices[-1]
        commandline = 'gmx_mpi grompp -f {} -c {} -p {} -n {} -o {} -maxwarn 3'.format(input_file, gro_file,
                                                                                             top_file, index_file,
                                                                                             output_file)
    elif thread.current_type == 'system_npt':
        input_file = get_input('npt_sys.mdp')
        index_file = thread.history.indices[-1]
        commandline = 'gmx_mpi grompp -f {} -c {} -p {} -n {} -o {} -maxwarn 3'.format(input_file, gro_file,
                                                                                             top_file, index_file,
                                                                                             output_file)
    elif thread.current_type == 'system_nvt':
        input_file = get_input('nvt_sys.mdp')
        index_file = thread.history.indices[-1]
        commandline = 'gmx_mpi grompp -f {} -c {} -p {} -n {} -o {} -maxwarn 3'.format(input_file, gro_file,
                                                                                             top_file, index_file,
                                                                                             output_file)
    else:
        pass

    subprocess.run(commandline, shell=True)
    # todo: try doing it this way and then do not need to have index file for peptide alone

def translate(thread, settings):
    # todo: testing needed
    gro_file = thread.history.coords[-1]
    output_file = thread.name + thread.suffix # todo: check!
    x, y = get_surface_size(settings) / 2 # todo: ensure both are divided by 2
    z = settings.initial_height
    commandline_arg = 'gmx_mpi editconf -f {} -o {} -translate {} {} {}'.format(gro_file, output_file, str(x), str(y), str(z))
    subprocess.run(commandline_arg, shell=True)
    thread.history.coords.append(output_file)

def convert_coord(thread, settings):
    # todo: testing needed
    # converting .pdb to .gro or .gro to .pdb
    gro_file = thread.history.coords[-1]
    pdb_file = thread.name # todo: check!
    commandline_arg = 'gmx_mpi editconf -f {} -o {}'.format(gro_file, pdb_file)
    subprocess.run(commandline_arg, shell=True)
    thread.history.coords.append(pdb_file)

def create_index(thread, settings):
    # todo: testing needed. Input file 'q' needed
    gro_file = thread.history.coords[-1]
    commandline = 'gmx_mpi make_ndx -f {} -o {}'.format(gro_file, index_file)
    subprocess.run(commandline, shell=True)
    thread.history.indices.append(index_file)

def combine_index(thread, settings):
    # todo: more pythonic way of doing this than cat!
    pass

def combine_pdb(thread):
    pass


