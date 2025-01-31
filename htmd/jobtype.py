"""
Interface for JobType objects. New JobTypes can be implemented by constructing a new class that inherits from JobType
and implements its abstract methods.
"""

import abc
import os
import sys
import time
import subprocess
import random
import pickle
import argparse
import numpy
import shutil
import time
#import pytraj
#import mdtraj
import warnings
import copy
import re
#import psutil
import utilities
import main
from infrastructure import factory


class JobType(abc.ABC):
    """
    Abstract base class for job types.
    Implements methods for all of the job type-specific tasks that ATESA might need.
    """

    @abc.abstractmethod
    def get_initial_coordinates(self, settings):
        """
        Obtain list of the appropriate initial coordinate files and copy them to the working directory.
        Parameters
        ----------
        settings : argparse.Namespace
            Settings namespace object
        Returns
        -------
        initial_coordinates : list
            List of strings naming the applicable initial coordinate files that were copied to the working directory
        """

        pass

    @abc.abstractmethod
    def get_next_step(self, thread, settings):
        """
        Determine and return name for next step in the thread given by "self"
        Parameters
        ----------
        thread : Thread()
            The Thread object on which to operate
        settings : argparse.Namespace
            Settings namespace object
        Returns
        -------
        name : str
            Name for the next step
        """

        pass

    @abc.abstractmethod
    def get_input_file(self, thread, settings):
        """
        Obtain appropriate input file for next job.
        At its most simple, implementations of this method can simply return settings.path_to_input_files + '/' +
        settings.job_type + '_' + settings.md_engine + '.in'
        Parameters
        ----------
        thread : Thread
            The Thread object on which to operate
        settings : argparse.Namespace
            Settings namespace object
        Returns
        -------
        input_file : str
            Name of the applicable input file
        """

        pass

    @abc.abstractmethod
    def get_batch_template(self, settings):
        """
        Return name of batch template file for the type of job indicated.
        Parameters
        ----------
        settings : argparse.Namespace
            Settings namespace object
        Returns
        -------
        name : str
            Name of the batch file template requested
        """

        pass

    @abc.abstractmethod
    def get_struct(self, thread):
        """
        Return the name of the appropriate inpcrd and topology files for the next step in the thread given by self.
        Parameters
        ----------
        thread : Thread
            The Thread object on which to operate
        Returns
        -------
        inpcrd : str
            String containing desired coordinate file name
        top : str
            String containing desired topology file name
        """

        pass

    @abc.abstractmethod
    def update_history(self, thread, settings, **kwargs):
        """
        Update or initialize the history namespace for this job type.
        This namespace is used to store the full history of a threads coordinate and trajectory files, as well as their
        results if necessary.
        If update_history is called with a kwargs containing {'initialize': True}, it simply prepares a blank
        history namespace and returns it. Otherwise, it adds the values of the desired keywords (which are desired
        depends on the implementation) to the corresponding history attributes in the index given by thread.suffix.
        Implementations of update_history should always save a copy of a freshly initialized thread.history object in a
        pickle file named 'algorithm_history.pkl' if one does not already exist. This is to support parallelization of
        algorithms across threads while still keeping the "shape" of the history object defined in one place (namely,
        in the implementation of update_history.)
        Parameters
        ----------
        thread : Thread
            The Thread object on which to operate
        settings : argparse.Namespace
            Settings namespace object
        kwargs : dict
            Dictionary of arguments that might be used to update the history object
        Returns
        -------
        None
        """

        pass

    @abc.abstractmethod
    def analyze(self, thread, settings):
        """
        Perform necessary analysis of a completed simulation step and store results as appropriate into thread.history
        Parameters
        ----------
        thread : Thread
            The Thread object on which to operate
        settings : argparse.Namespace
            Settings namespace object
        Returns
        -------
        None
        """

        pass

    @abc.abstractmethod
    def algorithm(self, thread, allthreads, settings):
        """
        Implement the algorithm that determines the next step and sets up the appropriate thread.history attributes
        Parameters
        ----------
        thread : Thread
            The Thread object on which to operate
        allthreads : list
            A list of all of the thread objects to consider during the algorithm.
        settings : argparse.Namespace
            Settings namespace object
        Returns
        -------
        terminate : Bool
            If True, terminate the entire isEE job
        """

        pass

    @abc.abstractmethod
    def gatekeeper(self, thread, allthreads, settings):
        """
        Return boolean indicating whether job is ready for next interpretation step.
        Parameters
        ----------
        thread : Thread
            The Thread object on which to operate
        allthreads : list
            List of all thread objects to consider
        settings : argparse.Namespace
            Settings namespace object
        Returns
        -------
        status : bool
            If True, ready for next interpretation step; otherwise, False
        """

        pass


class Adsorption_Simple(JobType):
    """
    Adapter class for primary Adsorption jobtype. My current conception is that this will be the only jobtype, but I've held
    off on simplifying the JobType class to only this implementation for now in case that changes
    """

    def get_initial_coordinates(self, settings):
        pass

    def get_next_step(self, thread, settings):
        # todo: will probably need to tweak this but a general idea
        if thread.current_type == '':  # this is the first step
            thread.current_type = 'peptide'
        elif thread.current_type == 'peptide':
            thread.current_type = 'system'
        elif thread.current_type == 'system':  # move on to next peptide
            if thread.current_peptide < len(thread.peptides)-1:  # should not ever fail this
                thread.current_type = 'peptide'
                thread.current_peptide += 1 # todo: what will it do when it gets to the end of the list - need it to stop running jobs for this thread
            else: # done processing peptides. send termination crierion
                thread.current_type = 'terminate'

        # todo: need to decide how it recognizes thread finished
        # todo: should this return something that will be equal to the thread.name in process.py?
        # todo: may be able to overlap some of this naming if essentially doing the same thing for nvt and npt peptide vs system
        return thread.current_type

    def get_input_file(self, thread, settings):
        return settings.path_to_input_files #+ '/' + settings.job_type + '_' + settings.md_engine + '.in'
        # todo: check this is applicable

    def get_batch_template(self, thread, settings):
        templ = settings.md_engine + '_' + settings.batch_system + '_' + thread.current_type + '.tpl'
        if os.path.exists(settings.path_to_templates + '/' + templ):
            return templ
        else:
            raise FileNotFoundError('cannot find required template file: ' + templ)
        # todo: check this is applicable

    def get_struct(self, thread):
        return thread.history.coords[thread.current_peptide][-1], thread.history.tops[thread.current_peptide][-1], \
               thread.history.indices[thread.current_peptide][-1]

    def update_history(self, thread, settings, **kwargs):
        if 'initialize' in kwargs.keys():
            if kwargs['initialize']: # todo: should these all actually be lists of lists for each peptide? - easier to extract data later...
                thread.history = argparse.Namespace() # todo: for now, just going to track current peptide w/ thread.peptide
                thread.history.peptides = []  # list of list of strings; initialized by main.init_threads(), updated by algorithm
                thread.history.trajs = []  # list of strings; updated by update_history() called by process.py (.xtc files)
                thread.history.tops = [[] for i in range(len(thread.peptides))]  # list of strings; initialized by main.init_threads(), updated by algorithm (.top filses)
                thread.history.coords = [[] for i in range(len(thread.peptides))]  # list of strings; initialized by main.init_threads(), updated by algorithm (.gro files)
                thread.history.indices = [[] for i in range(len(thread.peptides))]  # list of strings; initialized by main.init_threads(), updated by algorithm (.ndx files)
                thread.history.runfiles = [[] for i in range(len(thread.peptides))]  # list of strings; initialized by main.init_threads(), updated by algorithm (.tpr files)
                thread.history.ff = [[] for i in range(len(thread.peptides))]  # list of strings; initialized by main.init_threads(), updated by algorithm (.itp files)
                thread.history.timestamps = [[] for i in range(len(thread.peptides))]  # list of ints representing seconds since the epoch for the end of each step; initialized by main.init_threads(), updated by algorithm

            # todo: do I want to also keep track of .edr, .log, .tpr files?
            #if not os.path.exists( # todo: may not need timestamps?
            #        settings.working_directory + '/algorithm_history.pkl'):  # initialize algorithm_history file if necessary # todo: deprecate?
            #    pickle.dump(thread.history, open(settings.working_directory + '/algorithm_history.pkl',
            #                                     'wb'))  # an empty thread.history template

            # todo: this is probably not necessary. do need to update history.peptides for which have been processed
            #if 'add_peptides' in kwargs.keys():
            #    thread.history.peptides.extend(kwargs['add_peptides'])
        else:  # thread.history should already exist
            # todo: do I want to add trajectories in now? (npt.xtc, nvt.xtc created at each job step)

            # Add coordinate files from most recent batch job
            thread.history.coords[thread.current_peptide].append(kwargs['name'] + '_npt.gro')
            thread.history.coords[thread.current_peptide].append(kwargs['name'] + '_nvt.gro')

            # Add tpr files from the most recent batch job
            thread.history.runfiles[thread.current_peptide].append(kwargs['name'] + '_npt.tpr')
            thread.history.runfiles[thread.current_peptide].append(kwargs['name'] + '_nvt.tpr')

            # Add xtc files from the most recent batch job
            thread.history.runfiles[thread.current_peptide].append(kwargs['name'] + '_npt.xtc')
            thread.history.runfiles[thread.current_peptide].append(kwargs['name'] + '_nvt.xtc')


    def analyze(self, thread, settings):
        # If current_type == 'peptide', no analysis needed

        # If current_type == 'system', may need to do check for adsorption event
        pass

    def algorithm(self, thread, allthreads, running, settings):
        #this_algorithm = factory.algorithm_factory(settings.algorithm) # todo: not using algorithm class here...

        thread.current_type = thread.get_next_step(settings)
        # todo: implement: thread.current_peptide = thread.get_next_peptide(thread, settings)

        # If current system is only peptide
        if thread.current_type == 'peptide':

            # Build peptide in Amber TLEaP. This will add coord to history
            tleap_pdb = utilities.build_peptide(thread, settings) # todo: may or may not be returning...

            # Make pdb compatible with gromacs
            tleap_pdb_mod = utilities.edit_pdb(thread, settings)

            # Convert pdb to gro file
            peptide_gro = utilities.pdb2gmx(thread, settings)

            # Edit protein .itp file
            peptide_ff = utilities.clean_peptide_ff(thread, settings)

            # Edit topology template file to include peptide
            peptide_topology = utilities.write_topology(thread, settings)
            
            # Set box size and center peptide
            peptide_center = utilities.center_peptide(thread, settings)
            
            # Solvate box
            peptide_solvate = utilities.solvate(thread, settings)

            # Generate ion run file and ionize
            utilities.grompp_ion_runfile(thread, settings)
            peptide_ions = utilities.genion(thread, settings)

            # Add empty string to thread.history.index for peptide template
            thread.history.indices[thread.current_peptide].append('')

            # Ready to submit batch job!
            # todo: do I need to return anything? - will return a False termination criterion

        elif thread.current_type == 'system': # peptide and surface system

            # todo: for now, peptide starting coordinate will be nvt.gro (in future may use analyze to extract
            #  configuration)

            # Extract peptide from nvt trajectory
            peptide_relaxed = utilities.extract_peptide(thread, settings)

            # Grow peptide box and center at initial height
            system_center = utilities.center_peptide(thread, settings)

            # Edit topology template file to include peptide
            system_topology = utilities.write_topology(thread, settings)

            # Solvate box
            system_solvate = utilities.solvate(thread, settings)

            # Generate ion run file and ionize
            utilities.grompp_ion_runfile(thread, settings)
            system_ions = utilities.genion(thread, settings)

            # Translate box in z-direction to insert surface below
            system_translate = utilities.translate(thread, settings)

            # Convert peptide system coordinate file from .gro to .pdb
            peptide_gro_file = thread.history.coords[thread.current_peptide][-1]
            peptide_pdb_file = thread.name + '_' + thread.current_peptide + '_' + thread.current_type + '_trans.pdb'
            utilities.convert_coord(peptide_gro_file, peptide_pdb_file, thread, settings)
            # todo: need to see if easier in python to deal w/ pdb files or gro files (prob pdb)
            # todo: then probably want to make a surface .pdb only once - but need for entire thread, can store in thread.history

            # Combine system .pdb and surface .pdb
            surface_pdb = thread.surface_pbd
            combined_pdb = utilities.combine_pdb(surface_pdb, peptide_pdb_file, thread, settings)

            # Convert system coordinate file from .pdb to .gro
            combined_gro_file = thread.name + '_' + thread.current_peptide + '_comdbined.gro'
            utilities.convert_coord(combined_pdb, combined_gro_file, thread, settings)

            # Add surface ff to topology file
            combined_topology = utilities.add_to_topology(thread, settings)

            # Create index file for system and if slab center frozen combine to create new system index
            initial_index = utilities.create_index(thread, settings)
            if settings.frozen == True:
                system_index = utilities.combine_index(thread, settings)

            # Ready to submit batch job!

        # If the thread is done processing peptides, return a termination criterion for this thread only
        else:  # thread.current_type == 'terminate' (thread.current_peptide == len(thread.peptides) - 1):
            thread.terminate = True

        return False #running # todo: based on main loop, seems like this needs to return False if not terminating

    def gatekeeper(self, thread, allthreads, settings):
        # todo: consider adding a function w/in gatekeeper that can identify an adsorption event and bring about early termination - or would this be analysis?

        # If job for this thread has status 'C'ompleted/'C'anceled...
        if thread.get_status(0, settings) == 'C':  # index 0 because there is only ever one element in thread.jobids
            # todo: implement restarting if simulation crashed before a certain number of steps completed?
            return True
        else:
            return False

class Adsorption(JobType):
    """
    Adapter class for primary Adsorption jobtype. My current conception is that this will be the only jobtype, but I've held
    off on simplifying the JobType class to only this implementation for now in case that changes
    """

    def get_initial_coordinates(self, settings):
        pass

    def get_next_step(self, thread, settings):
        # Number of steps for "system" will depend on inputted number of structures and replicates
        if thread.current_type == '':  # this is the first step
            thread.current_type = 'peptide'
        elif thread.current_type == 'peptide':
            thread.current_type = 'system'
        elif thread.current_type == 'system':
            # If system is built then move on to replicates with same system
            thread.current_type = 'replicate'
        elif thread.current_type == 'replicate':
            # Repeat system run for number of replicates
            if int(thread.current_rep) < settings.num_reps - 1:
                thread.current_rep += 1
            elif int(thread.current_struct) < settings.num_structures - 1: # finished the last replicate for the previous structure
                thread.current_struct += 1
                thread.current_rep = 0
                thread.current_type = 'system'
            elif int(thread.current_peptide) < len(thread.peptides) - 1: # finished the last structure for the previous peptide
                thread.current_struct = 0
                thread.current_rep = 0
                thread.current_type = 'peptide'
            else:
                thread.current_type = 'terminate'
        return thread.current_type


        #elif thread.current_type == 'system':
            # starting with thread.current_struct == 0, thread.current_rep == 0:
        #    if int(thread.current_rep) < settings.num_reps - 1:
        #        thread.current_rep += 1
        #    elif int(thread.current_struct) < settings.num_structures - 1: # just completed the last replicate for the previous structure
        #        thread.current_struct += 1
        #    elif thread.current_peptide < len(thread.peptides) - 1: # finished the last structure for the previous peptide
        #        thread.current_peptide += 1
        #        thread.current_type = 'peptide'
        #    else:  # done processing peptides. Send termination criterion
        #        thread.current_type = 'terminate' # todo: is this the correct way to say that thread is done being processed?

    def get_input_file(self, thread, settings):
        return settings.path_to_input_files #+ '/' + settings.job_type + '_' + settings.md_engine + '.in'
        # todo: check this is applicable

    def get_batch_template(self, thread, settings):
        templ = settings.md_engine + '_' + settings.batch_system + '_' + thread.current_type + '.tpl'
        if os.path.exists(settings.path_to_templates + '/' + templ):
            return templ
        else:
            raise FileNotFoundError('cannot find required template file: ' + templ)
        # todo: check this is applicable

    def get_struct(self, thread):
        if thread.current_type == 'peptide':
            return thread.history.pep_coords[thread.current_peptide], thread.history.pep_tops[thread.current_peptide], \
               thread.history.pep_indices[thread.current_peptide]
        else: # thread.current_type == 'system' or 'replicate'
            return thread.history.coords[thread.current_peptide][thread.current_struct], thread.history.tops[thread.current_peptide][thread.current_struct], \
                thread.history.indices[thread.current_peptide][thread.current_struct]


    def update_history(self, thread, settings, **kwargs):
        if 'initialize' in kwargs.keys():
            if kwargs['initialize']: # todo: should these all actually be lists of lists for each peptide? - easier to extract data later...

                # Peptide history attributes handle variable number of peptides
                thread.history = argparse.Namespace()  # todo: for now, just going to track current peptide w/ thread.peptide
                thread.history.peptides = []  # list of list of strings; initialized by main.init_threads(), updated by algorithm
                thread.history.pep_trajs = [] # list of strings where index corresponds to current peptide
                thread.history.pep_tops = []
                thread.history.pep_coords = []
                thread.history.pep_indices = []
                thread.history.pep_npt = []
                thread.history.pep_nvt = []
                thread.history.pep_ff = []

                # populated at end of peptide run
                thread.history.final_coords = [] # todo: keep track of these or just go off of file naming convention?
                thread.history.final_xtc = []
                thread.history.final_tpr = []

                # populated during analysis
                thread.history.structs = [] # todo: make sure this works - will be list of lists eventually

                # System history attributes handle variable number of peptides, structures, and replicates
                thread.history.tops = [[] for i in range(len(thread.peptides))]
                thread.history.coords = [[] for i in range(len(thread.peptides))]
                thread.history.indices = [[] for i in range(len(thread.peptides))]

                # thread.history.trajs = [[[[] for k in range(settings.num_reps)] for j in range(settings.num_structures)] for i in range(len(thread.peptides))]  # (.xtc files)
                #thread.history.tops = [[[[] for k in range(settings.num_reps)] for j in range(settings.num_structures)] for i in range(len(thread.peptides))]  # (.top files)
                #thread.history.coords = [[[[] for k in range(settings.num_reps)] for j in range(settings.num_structures)] for i in range(len(thread.peptides))]  # (.gro files)
                #thread.history.indices = [[[[] for k in range(settings.num_reps)] for j in range(settings.num_structures)] for i in range(len(thread.peptides))]  # (.ndx files)
                #thread.history.runfiles = [[[[] for k in range(settings.num_reps)] for j in range(settings.num_structures)] for i in range(len(thread.peptides))]  # (.tpr files)
                #thread.history.ff = [[[[] for k in range(settings.num_reps)] for j in range(settings.num_structures)] for i in range(len(thread.peptides))]  # (.itp files)
                # list of ints representing seconds since the epoch for the end of each step; initialized by main.init_threads(), updated by algorithm

        else:  # thread.history should already exist
            # todo: need to check, but peptide should just be acting on struct,rep = 0

            # todo: need to make consistent with history object naming based on previous jobtype

            if thread.current_type == 'peptide':

                # Add coordinate files from most recent batch job
                # todo: need to make consistent with history object naming based on previous jobtype

                # Add coordinate files from most recent batch job
                thread.history.final_coords.append(kwargs['name'] + '_npt.gro')
                #thread.history.coords[thread.current_peptide][thread.current_struct][thread.current_rep].append(kwargs['name'] + '_npt.gro')
                #thread.history.coords[thread.current_peptide][thread.current_struct][thread.current_rep].append(kwargs['name'] + '_nvt.gro')

                # Add tpr files from the most recent batch job
                thread.history.final_xtc.append(kwargs['name'] + '_nvt.tpr')
                #thread.history.runfiles[thread.current_peptide][thread.current_struct][thread.current_rep].append(kwargs['name'] + '_npt.tpr')
                #thread.history.runfiles[thread.current_peptide][thread.current_struct][thread.current_rep].append(kwargs['name'] + '_nvt.tpr')

                # Add xtc files from the most recent batch job
                thread.history.final_tpr.append(kwargs['name'] + '_nvt.xtc')
                #thread.history.trajs[thread.current_peptide][thread.current_struct][thread.current_rep].append(kwargs['name'] + '_npt.xtc')
                #thread.history.trajs[thread.current_peptide][thread.current_struct][thread.current_rep].append(kwargs['name'] + '_nvt.xtc')

            else:
                pass # todo: nothing being updated at end of other runs?




    def analyze(self, thread, settings):
        # If current_type == 'peptide' --> pulling out structures here

        # todo: add structure to thread.history.coords[thread.current_peptide][thread.current_stuct][for all thread.current_rep]
        # todo: that way all reps start with the same initial peptide structure

        if thread.current_type == 'peptide':

            # randomly sample from nvt trajectory to create intial coordinate files for num_structures, common among num_reps
            coord = thread.history.final_coords[thread.current_peptide]
            xtc = thread.history.final_xtc[thread.current_peptide]
            tpr = thread.history.final_tpr[thread.current_peptide]
            coord_files = utilities.sample_trajectory(thread, settings, xtc, tpr)

            # add coordinate files to history object
            thread.history.structs.append(coord_files)

        else: # thread.current_type == system
            pass

    def algorithm(self, thread, allthreads, running, settings):

        thread.current_type = thread.get_next_step(settings)
        # todo: implement: thread.current_peptide = thread.get_next_peptide(thread, settings)

        # If current system is peptide only
        if thread.current_type == 'peptide':
            print('beginning peptide build')
            # Build peptide in Amber TLEaP. This will add coord to history
            tleap_pdb = utilities.build_peptide(thread, settings)  # todo: may or may not be returning...
            print('built tleap pdb')
            # Make pdb compatible with gromacs
            tleap_pdb_mod = utilities.edit_pdb(thread, settings, tleap_pdb)
            print('pbd made gmx compatible')
            # Convert pdb to gro file
            peptide_gro = utilities.pdb2gmx(thread, settings, tleap_pdb_mod)
            print('pdb converted to gmx')
            # Edit protein .itp file
            peptide_ff_in = thread.name + '_' + str(thread.current_peptide) + '.top'
            print(os.path.exists(settings.working_directory + peptide_ff_in))
            print(peptide_ff_in)
            peptide_ff = utilities.clean_peptide_ff(thread, settings, peptide_ff_in)

            # Edit topology template file to include peptide
            peptide_topology = utilities.write_topology(thread, settings, peptide_ff)

            # Set box size and center peptide
            peptide_center_gro = utilities.center_peptide(thread, settings, peptide_gro)

            # Solvate box
            peptide_solvate_gro = utilities.solvate(thread, settings, peptide_center_gro, peptide_topology)

            # Generate ion run file and ionize
            tpr_file = utilities.grompp_ion_runfile(thread, settings, peptide_solvate_gro, peptide_topology)
            peptide_ion_gro = utilities.genion(thread, settings, tpr_file, peptide_topology, peptide_ff)

            # Add empty string to thread.history.index for peptide template
            thread.history.pep_indices.append('') # todo: actually I don't think I need this anymore?

            # Ready to submit batch job! Add simulation ready files to thread.history namespace
            # Needed - final .gro file, .top file, .ndx file
            thread.history.pep_coords.append(peptide_ion_gro)
            thread.history.pep_tops.append(peptide_topology)
            # todo: do I need to return anything? - will return a False termination criterion

        elif thread.current_type == 'system': # peptide and surface system

            print('SYSTEM STEP')

            # Grow peptide box and place at randomly selected x-y coordinates within the range of the surface
            peptide_coord = thread.history.structs[thread.current_peptide][thread.current_struct]
            place_peptide_gro = utilities.place_peptide(thread, settings, peptide_coord)
            # todo: this step is done for each rep, but would actually be doing the same thing for each rep (could move this to analysis?)

            print(thread.history.pep_ff, thread.current_peptide)
            # Edit topology template file to include peptide
            peptide_ff = thread.history.pep_ff[thread.current_peptide]
            system_topology = utilities.write_topology(thread, settings, peptide_ff)

            # Solvate box
            system_solvate_gro = utilities.solvate(thread, settings, place_peptide_gro, system_topology)

            # Generate ion run file and ionize
            tpr_file = utilities.grompp_ion_runfile(thread, settings, system_solvate_gro, system_topology)
            system_ions_gro = utilities.genion(thread, settings, tpr_file, system_topology)

            # Translate box in z-direction to insert surface below
            system_translate_gro = utilities.translate(thread, settings, system_ions_gro)

            # Convert peptide system coordinate file from .gro to .pdb
            peptide_gro_file = thread.history.coords[thread.current_peptide][-1]
            peptide_pdb_file = thread.name + '_' + thread.current_peptide + '_' + thread.current_type + '_trans.pdb'
            peptide_pdb = utilities.convert_coord(thread, settings, peptide_gro_file, peptide_pdb_file)
            # todo: need to see if easier in python to deal w/ pdb files or gro files (prob pdb)
            # todo: then probably want to make a surface .pdb only once - but need for entire thread, can store in thread.history

            # Combine system .pdb and surface .pdb
            surface_pdb = thread.surface_pbd
            combined_pdb = utilities.combine_pdb(thread, settings, peptide_pdb, surface_pdb)

            # Convert system coordinate file from .pdb to .gro
            combined_gro_file = thread.name + '_' + thread.current_peptide + '_comdbined.gro'
            utilities.convert_coord(thread, settings, combined_pdb, combined_gro_file)

            # Add surface ff to topology file
            combined_topology = utilities.add_to_topology(thread, settings, system_topology)

            # Create index file for system and if slab center frozen combine to create new system index
            system_index = utilities.create_index(thread, settings, combined_gro_file)
            if settings.frozen == True:
                system_index = utilities.combine_index(thread, settings, system_index)

            # Ready to submit batch job!
            # Add simulation ready files to thread.history namespace
            thread.history.tops.append(system_topology)
            thread.history.coords.append(combined_gro_file)
            thread.history.indices.append(system_index)


        elif thread.current_type == 'replicate':

            # using topology, coordinate, and index files for the current structure and peptide
            pass

        else: # thread.current_type == 'terminate'

            # todo: do I need to set a termination criteria here?
            pass



        # Intermediate steps:

        # 1. Build peptide in water system (same as Simple Adsorption_Simple jobtype
        # current_type == "peptide"

        # 2. Energy minimize and equilibrate peptide system

        # 3. Extract 4 structures from nvt trajectory of peptide
        # Done within the analyze step

        # 4. Graft peptide onto slab
        # current_type == "system"
        # done for each of the 4 structures
        # grafting onto 10 different random locations of the surface

        # 5. Energy minimize and equilibrate peptide-surface systems


    def gatekeeper(self, thread, allthreads, settings):
        # todo: consider adding a function w/in gatekeeper that can identify an adsorption event and bring about early termination - or would this be analysis?

        # If job for this thread has status 'C'ompleted/'C'anceled...
        if thread.get_status(0, settings) == 'C':  # index 0 because there is only ever one element in thread.jobids
            # todo: implement restarting if simulation crashed before a certain number of steps completed?
            return True
        else:
            return False
