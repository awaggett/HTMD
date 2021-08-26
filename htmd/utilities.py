"""
Utility functions implemented here are clearly defined unit operations. They may only be called once in the code, but
are defined separately for cleanliness and legibility.
"""

import os
import re
import shutil
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

def get_input(settings, input_file):
    if os.path.exists(settings.path_to_inout_files + '/' + input_file):
        return os.path.join(settings.path_to_inout_files, input_file)  # todo: do I need the full path here? or return input_file?

def get_template(settings, templ):
    if os.path.exists(settings.path_to_templates + '/' + templ):
        return os.path.join(settings.path_to_templates, templ)  # todo: do I need the full path here?

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
    # write sequence in tleap format, e.g. ['ACE', 'LYS', 'NME'] --> '{ ACE LYS NME }
    sequence = '{ ' + " ".join(thread.peptides[thread.current_peptide]) + ' }'
    # or for bash script purposes: input seq = thread.peptides[thread.peptide]
    # sequence = '{ ' + " ".join(seq) + ' }'

    try:
        system = tleap()
    except TypeError:
        system = tleap.System()
    short_name = thread.name + '_' + thread.current_peptide
    long_name = short_name + '_tleap.pdb'
    system.pbc_type = None
    system.neutralize = False
    system.output_path = settings.working_directory # todo: check on this
    system.output_prefix = thread.name + '_' + thread.current_peptide + '_tleap' # todo: will need to specify which peptide also
    system.template_lines = [
        'source oldff/leaprc.ff99SBildn', # todo: update ff, not completely necessary
        short_name + ' = sequence ' + sequence, # todo: sequence = thread.current_peptides[i]
        'savepdb ' + short_name + ' ' + long_name, # todo: decide on file naming
        'quit',
    ]

    # add output coordinate file to thread.history.coords
    thread.history.coords[thread.current_peptide].append(long_name)
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
    pdb_in = thread.history.coords[thread.current_peptide][-1]
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
            aa = line_split[6]
            res = line_split[8]

            # change histidine naming
            if aa == 'HIE' or aa == 'HIP' or aa == 'HID': # todo: determine how best to deal w/ protonation
                line_split[6] = 'HIS'

            # special case for  SER, CYS
            if (aa == 'CYS' or aa == 'SER') and atom == 'HG':
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
        thread.history.coords[thread.current_peptide].append(pdb_out)
    return pdb_out

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
    pdb_file = thread.history.coords[thread.current_peptide][-1]
    gro_file = thread.name + '_' + thread.current_peptide + '_init.gro'
    protein_ff = settings.force_field.split('.')[0]
    topol_file = thread.name + '_' + thread.current_peptide + '.top'
    posre_file = thread.name + '_' + thread.current_peptide + '_posre.itp'

    commandline_arg = 'echo 3 4 | gmx_mpi pdb2gmx -f {} -o {} -ff {} -p {} -i {} -water none -ter '.format(pdb_file, gro_file, protein_ff, topol_file, posre_file)
    subprocess.run(commandline_arg, shell=True)

    thread.history.coords[thread.current_peptide].append(gro_file)
    thread.history.peptide_ff[thread.current_peptide].append(topol_file)  # topology file will be edited and replaced with .itp file

    # remove posre.itp - will not be used # todo: in the future may want to allow for option to keep
    #commandline_arg2 = 'rm {}'.format(posre_file)
    #subprocess.run(commandline_arg2, shell=True)
    os.remove(posre_file)

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

    protein_ff_in = thread.history.tops[thread.current_peptide][-1]
    protein_ff_out = thread.name + '_' + thread.current_peptide + '.itp'
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
            elif start_writing is True and end_writing is False:
                ff_out.write(line)

    # add peptide.itp to history namespace
    thread.history.peptide_ff[thread.current_peptide].append(protein_ff_out)

    # remove gromacs generated topology file and posre file
    #os.remove('posre.itp') # tood: done in pdb2gmx, but could do here instead?
    os.replace(protein_ff_in, protein_ff_out) # todo: better to just remove?

    return protein_ff_out

def write_topology(thread, settings):
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
    ff = settings.force_field
    peptide_ff = thread.history.ff[thread.current_peptide][-1]
    water_ff = settings.path_to_water_ff
    ion_ff = settings.path_to_ion_ff
    thread_name = thread.name
    thread_peptide = thread.current_peptide

    lines = ["; Include forcefield parameters\n", "#include \"./" + str(ff) + "\"\n", "#include \"./" + str(peptide_ff)
             + "\"\n", "\n", "; Include water topology\n", "#include \"./" + str(water_ff) + "\"\n", "\n",
             "; Include topology for ions\n", "#include \"./" + str(ion_ff) + "\"\n", "\n", "[ system ]\n", "; Name\n",
             str(thread_name) + " " + str(thread_peptide) + "\n", "\n", "[ molecules ]\n", "; Compound        #mols\n",
             "\n", "Protein              1"]

    topol_out = open(thread.name + '_' + thread.current_peptide + '_' + thread.current_type + '_topol.top', "a")
    topol_out.writelines(lines)
    topol_out.close
    thread.history.tops.[thread.current_peptide].append(topol_out)
    # add protein ff to topology file
    #with open(topol_out, 'rw'):
    #    for line in topol_out:
    #        line_split = re.split(r'(\s+)', line)
    #        if line_split[0] == '#protein_ff':
    #            top_out.write('#include "./{}.itp"'.format(base_name))
    #        elif line_split[0] == 'sys_name':
    #            top_out.write('{} on Quartz slab'.format(base_name))
    #        elif line_split[0] == 'protein_name':
    #            top_out.write('{}       1'.format(base_name))
    #        else:
    #            top_out.write(line)

    return topol_out

def get_surface_size(settings):
    surface = settings.surface_coord
    dims = open(surface, 'r').readlines()[-1].split()
    return dims[0], dims[1], dims[2]
    # todo: check this returns the correct thing - could also just pull from .itp file

def center_peptide(thread, settings):
    # todo: testing needed !!!
    gro_file = thread.history.coords[thread.current_peptide][-1]
    output_file = thread.name + '_' + thread.current_peptide + '_' + thread.current_type + '_center.gro' # todo: make sure these are all strings...

    if thread.current_type == 'peptide':
        x = settings.peptide_box_x
        y = settings.peptide_box_y
        z = settings.peptide_box_z
        commandline_arg = 'gmx_mpi editconf -f {} -o {} -box {} {} {} '.format(gro_file, output_file, x, y, z)
    else:  # thread.current_type == 'system'
        x, y, z = get_surface_size(settings)
        z_mod = z - settings.slab_height
        peptide_x = x/2 # todo: this is not compatible of non-cubic boxes
        peptide_y = y/2
        peptide_z = settings.initial_height
        commandline_arg = 'gmx_mpi editconf -f {} -o {} -box {} {} {} -center {} {} {} '.format(gro_file, output_file, x, y, z_mod, peptide_x, peptide_y, peptide_z)

    #gmx_mpi editconf -f {} -o {} -box {} {} {} '.format(gro_file, output_file, x, y, z)
    subprocess.run(commandline_arg, shell=True)
    thread.history.coords[thread.current_peptide].append(output_file)
    return output_file

def solvate(thread, settings):
    # todo: testing needed
    gro_file = thread.history.coords[thread.current_peptide][-1]
    output_file = thread.name + '_' + thread.current_peptide + '_' + thread.current_type + '_solvate.gro'
    topol_file = thread.history.tops[thread.current_peptide][-1]
    commandline_arg = 'gmx_mpi solvate -cp {} -cs spc216.gro -o {} -p {}'.format(gro_file, output_file, topol_file)
    subprocess.run(commandline_arg, shell=True)
    thread.history.coords[thread.current_peptide].append(output_file) # topology file name will be unchanged

def grompp_ion_runfile(thread, settings):
    # todo testing needed
    gro_file = thread.history.coords[thread.current_peptide][-1]
    tpr_file = thread.name + '_' + thread.current_peptide + '_' + thread.current_type + '_ion.tpr'
    topol_file = thread.history.tops[thread.current_peptide][-1]
    commandline_arg = 'gmx_mpi grompp -f ion.mdp -c {} -p {} -o {} -maxwarn 4'.format(gro_file, topol_file, tpr_file)
    subprocess.run(commandline_arg, shell=True)
    thread.history.runfiles[thread.current_peptide].append(tpr_file)

def get_system_charge(thread):
    # read charge from peptide .itp file
    itp_file = thread.history.ff[thread.current_peptide][-1]

    with open(itp_file, 'r') as itp: # todo: make sure this is always the case
        final_charge = int([line for line in itp.readlines() if 'qtot' in line][-1].split()[-1])
    return final_charge

def genion(thread, settings):
    # todo testing needed
    tpr_file = thread.history.runfiles[thread.current_peptide][-1]
    output_file = thread.name + '_' + thread.current_peptide + '_' + thread.current_type + '_ion.gro' # todo: check!
    charge = get_system_charge(thread)
    if charge > 0:
        charge_mod = '-np 5'
    else:  # charge <= 0
        charge_mod = '-nn 5' # todo: check if this works as expected for system
    commandline_arg = 'echo SOL | gmx_mpi genion -s {} -o {} -p topol.top -pname NA -nname CL -neutral {}'.format(tpr_file, output_file, charge_mod)
    subprocess.run(commandline_arg, shell=True)
    thread.history.coords[thread.current_peptide].append(output_file)
    return output_file

def extract_peptide(thread, settings):
    gro_file = thread.history.coords[thread.current_peptide][-1]
    tpr_file = thread.history.runfiles[thread.current_peptide][-1]
    output_file = thread.name + '_' + thread.current_peptide + '_peptide_relax.gro'

    # todo: will probably need to do this w/ POpen to read lines
    commandline_arg = 'echo 1 | gmx_mpi trjconv -f {} -s {} -o {} -pbc whole'.format(gro_file, tpr_file, output_file)
    subprocess.run(commandline_arg, shell=True)
    thread.history.coords[thread.current_peptide].append(output_file)

def translate(thread, settings):
    # todo: testing needed
    gro_file = thread.history.coords[thread.current_peptide][-1]
    output_file = thread.name + '_' + thread.current_peptide + '_' + thread.current_type + '_trans.gro'
    height = settings.initial_height

    commandline_arg = 'gmx_mpi editconf -f {} -o {} -translate 0 0 {}'.format(gro_file, output_file, height)
    subprocess.run(commandline_arg, shell=True)
    thread.history.coords[thread.current_peptide].append(output_file)

def combine_pdb(surface, peptide, thread, settings):
    output_file = thread.name + '_' + thread.current_peptide + '_comdbined.pdb'

    fout = open(output_file, 'w')
    fout.writelines(open(surface, 'r').readlines()[:-2])
    fout.writelines(open(peptide, 'r').readlines()[4:])

    thread.history.coords[thread.current_peptide].append(output_file)
    return output_file

def convert_coord(file_in, file_out, thread, settings):
    commandline_arg = 'gmx_mpi editconf -f {} -o {}'.format(file_in, file_out)
    subprocess.run(commandline_arg, shell=True)
    thread.history.coords[thread.current_peptide].append(file_out)

def add_to_topology(thread, settings):
    topology_file = thread.history.tops[thread.current_peptide][-1]
    surface_name = settings.surface_name
    output_file = thread.name + '_' + thread.current_peptide + '_' + thread.current_type + '_topol_sys.top'
    surface_topology = settings.surface_ff

    fout = open(output_file, 'w')
    lines = open(topology_file, 'r').readlines()
    for i, line in enumerate(lines):
        fout.write(line)
        if i == 1:
            fout.write('#include \"./' + surface_topology + '\"\n')
        if i == 15:
            fout.write(surface_name + '\t\t1\n')

    thread.history.tops[thread.current_peptide].append(output_file)
    return output_file

def create_index(thread, settings):
    # todo: testing needed. Input file 'q' needed
    gro_file = thread.history.coords[thread.current_peptide][-1]
    index_file = thread.name + '_' + thread.current_peptide + '_' + thread.current_type + "_init.ndx"

    commandline = 'echo q | gmx_mpi make_ndx -f {} -o {}'.format(gro_file, index_file)
    subprocess.run(commandline, shell=True)
    thread.history.indices[thread.current_peptide].append(index_file)
    return index_file

def combine_index(thread, settings):
    # todo: more pythonic way of doing this than cat!
    init_index = thread.history.indices[thread.current_peptide][-1]
    slab_center = settings.center
    output_file = thread.name + '_' + thread.current_peptide + '_' + thread.current_type + "_system.ndx"

    output = open(output_file, 'w')
    output.write(open(init_index, 'r').read())
    output.write(open(slab_center, 'r').read())

    thread.history.indices[thread.current_peptide].append(output_file)
    return output_file

def combine_pdb(thread):
    pass


