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
#try:
#    from paprika.build.system import TLeap as tleap
#except ModuleNotFoundError:
#    from paprika import tleap



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
    init_coords or None?

    """
    # todo: Call to TLEaP. Determine where/how files will be outputted. Either return file name of .pdb file generated
    # todo: or directly add to settings.init_coord(?) w/in this function

    pass

def edit_pdb(thread, settings ):
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
    pass

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
    pass

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
    pass

def center_peptide(thread, settings):
    pass


