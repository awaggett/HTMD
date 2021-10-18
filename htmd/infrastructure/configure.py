"""
configure.py
Takes user input file and returns settings namespace object
"""

import argparse
import pytraj       # to support pytraj calls in input file
import mdtraj       # to support mdtraj calls in the input file
import numpy        # to support numpy  calls in input file
import numpy as np  # to support numpy  calls even if called as np
import sys
import os
import shutil
import pickle
from jinja2 import Environment, FileSystemLoader
import typing
import pydantic

def configure(input_file, user_working_directory=''):
    """
    Configure the settings namespace based on the config file.

    Parameters
    ----------
    input_file : str
        Name of the configuration file to read
    user_working_directory : str
        User override for working directory (overrides value in input_file), ignored if set to ''

    Returns
    -------
    settings : argparse.Namespace
        Settings namespace object

    """

    class Settings(pydantic.BaseModel):
        # This class initializes the settings object with type hints. After being built, it gets exported as an
        # argparse.Namelist object, just for convenience.

        # Core settings required for all jobs
        job_type: str = 'adsorption'    # 'adsorption is currently the only jobtype written
        batch_system: str
        restart: bool           # todo: figure out unpickling
        md_engine: str = 'gromacs'      # 'gromacs' is currently the only MD engine adapted
        task_manager: str = 'simple'    # 'simple' is the only implementation I wrote
        working_directory: str
        overwrite: bool
        account: str 
        partition: str

        # Batch template settings
        nodes: int = 1                  # number of nodes will determine partitioning of threads
        ppn: int = 1
        mem: str = '20GB'
        walltime: str = '02:00:00'
        extra: str = ''

        # File path settings
        print(os.path.realpath(__file__))
        print(os.path.dirname(os.path.realpath(__file__)))
        path_to_input_files: str = os.path.dirname(os.path.realpath(__file__)) + '/data/input_files'
        path_to_templates: str = os.path.dirname(os.path.realpath(
            __file__)) + '/data/templates'  # todo: note that this is used below for establishing the Jinja2 environment


        # todo: replace the below and add your own as needed. The below lines are left in place only as an example.
        # Settings for Adsorption jobs
        name: str                       # identify job names
        peptide_box_x: int = 6          # x dimension of box to relax peptide in nm (default 6x6x6)
        peptide_box_y: int = 6          # y dimension of box
        peptide_box_z: int = 6          # z dimension of box
        force_field: str = 'charmm36-nov2016-repart.ff' # todo gromacs input split on .
        path_to_ion_ff: str = 'charmm36-nov2016-repart.ff/ions.itp'
        path_to_water_ff: str = 'charmm36-nov2016-repart.ff/spce_part.itp'
        initial_height: int = 2         # initial height of peptide above slab (default 2 nm)
        peptides: typing.List[typing.List[str]] = []
        surface_ff: str                 # .itp force field file for surface
        surface_coord: str              # coordinate file for surface
        surface_height: int             # approximate height (nm) of surface. Necessary for building peptide-surface system with sufficient space
        surface_pdb: str                # .pdb file for surface # todo: may want to have option to convert from .gro
        surface_name: str               # name of surface in .itp file
        surface_frozen: bool            # Boolean: True if center of slab is frozen
        surface_center: str             # txt file containing atom number of frozen center (needed if surface_frozen == True)
        # Batch template settings energy minimization (peptide or system)


    # Import config file line-by-line using exec()
    lines = open(input_file, 'r').readlines()
    line_index = 0
    for line in lines:      # each line in the input file is just python code setting a variable...
        line_index += 1
        try:
            exec(line)      # ... this means that comments are supported using '#' and whitespace is ignored.
        except Exception as e:
            raise ValueError('error raised while reading line ' + str(int(line_index)) + ' of configuration file '
                             + input_file + ': ' + str(e))

    # Define settings namespace to store all these variables
    config_dict = {}
    config_dict.update(locals())
    settings = argparse.Namespace()
    settings.__dict__.update(Settings(**config_dict))

    # Override working directory if provided with user_working_directory
    if user_working_directory:
        settings.working_directory = user_working_directory

    # Format directories properly (no trailing '/')
    if settings.working_directory[-1] == '/':
        settings.working_directory = settings.working_directory[:-1]
    if settings.path_to_input_files[-1] == '/':
        settings.path_to_input_files = settings.path_to_input_files[:-1]
    if settings.path_to_templates[-1] == '/':
        settings.path_to_templates = settings.path_to_templates[:-1]

    # Set Jinja2 template environment
    if os.path.exists(settings.path_to_templates):
        settings.env = Environment(loader=FileSystemLoader(settings.path_to_templates))
    else:
        raise FileNotFoundError('could not locate templates folder: ' + settings.path_to_templates)

    return settings

if __name__ == "__main__":
    configure('','')
