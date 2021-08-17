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
from htmd import utilities
from htmd import main
from htmd.infrastructure import factory


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


class Adsorption(JobType):
    """
    Adapter class for primary Adsorption jobtype. My current conception is that this will be the only jobtype, but I've held
    off on simplifying the JobType class to only this implementation for now in case that changes
    """

    def get_initial_coordinates(self, settings):
        pass

    def get_next_step(self, thread, settings):
        # todo: will probably need to tweak this but a general idea
        if thread.current_type == '': # this is the first step
            thread.current_type = 'peptide'
        elif thread.current_type == 'peptide':
            thread.current_type = 'system'
        elif thread.current_type == 'system': # move onto next peptide
            thread.current_type = 'peptide'
        # todo: need to decide how it recognizes thread finished
        # todo: should this return something that will be equal to the thread.name in process.py?
        # todo: may be able to overlap some of this naming if essentially doing the same thing for nvt and npt peptide vs system
        return thread.current_type

    def get_input_file(self, thread, settings):
        return settings.path_to_input_files + '/' + settings.job_type + '_' + settings.md_engine + '.in'
        # todo: check this is applicable

    def get_batch_template(self, settings):
        # todo: may change this template definition back.../ need to add more templates for system
        templ = settings.md_engine + '_' + settings.batch_system + '_' + thread.system + '_' + thread.current_type + '.tpl'
        if os.path.exists(settings.path_to_templates + '/' + templ):
            return templ
        else:
            raise FileNotFoundError('cannot find required template file: ' + templ)
        # todo: check this is applicable

    def get_struct(self, thread):
        #if settings.batch_system == 'gromacs':
        #    return thread.history.runfiles[-1] # todo: actually will be grompping w/in
        #else:
        return thread.history.coords[-1], thread.history.tops[-1], thread.history.index[-1]
        # todo: check this is applicable - may only need run file

    def update_history(self, thread, settings, **kwargs):
        if 'initialize' in kwargs.keys():
            if kwargs['initialize']:
                thread.history = argparse.Namespace()
                thread.history.peptides = []  # list of list of strings; initialized by main.init_threads(), updated by algorithm
                #thread.history.trajs = []  # list of strings; updated by update_history() called by process.py
                thread.history.tops = []  # list of strings; initialized by main.init_threads(), updated by algorithm
                thread.history.coords = []  # list of strings; initialized by main.init_threads(), updated by algorithm
                thread.history.indices = []  # list of strings; initialized by main.init_threads(), updated by algorithm
                thread.history.runfiles = []  # list of strings; initialized by main.init_threads(), updated by algorithm
                thread.history.ff = []  # list of strings; initialized by main.init_threads(), updated by algorithm
                thread.history.timestamps = []  # list of ints representing seconds since the epoch for the end of each step; initialized by main.init_threads(), updated by algorithm
            #if not os.path.exists( # todo: may not need timestamps?
            #        settings.working_directory + '/algorithm_history.pkl'):  # initialize algorithm_history file if necessary # todo: deprecate?
            #    pickle.dump(thread.history, open(settings.working_directory + '/algorithm_history.pkl',
            #                                     'wb'))  # an empty thread.history template

            # todo: this is probably not necessary. do need to update history.peptides for which have been processed
            #if 'add_peptides' in kwargs.keys():
            #    thread.history.peptides.extend(kwargs['add_peptides'])
        else:  # thread.history should already exist
            pass
            # todo: will probably need to put something here that removes peptides from thread.history once first process step is complete


    def analyze(self, thread, settings):
        pass

    def algorithm(self, thread, allthreads, running, settings):
        this_algorithm = factory.algorithm_factory(settings.algorithm)

        # todo: write if statement based on jobtype.current (=petide or system) may need to have current_type == '' then next step is peptide?
#        if thread.current_type == '':  # if this is the first step in this thread
#            next_step = this_algorithm.get_first_step(thread, allthreads, settings) # todo: may not need a return. can just update type
#        else:
#            next_step = this_algorithm.get_next_step(thread, allthreads, settings)
        thread.current_type = thread.get_next_step(thread, settings)

        # if current system is only peptide
        if thread.current_type == 'peptide':

            # build peptide in Amber TLEaP. This will add coord to history
            tleap_pdb = utilities.build_peptide(thread, settings) # todo: may or may not be returning...

            # make pdb compatible with gromacs
            tleap_pdb_mod = utilities.edit_pdb(thread, settings)

            # convert pdb to gro file
            peptide_gro = utilities.pdb2gmx(thread, settings)

            # edit protein .itp file and populate topology base file
            # todo: how to do this systematically

            # set box size and center peptide
            utilities.center_peptide(thread, settings)

            # solvate box
            utilities.solvate(thread, settings)

        else: # thread.current_type == 'system' (peptide and surface system)

            # add nvt.gro to history.coords?
            # center peptide

            # solvate box

            # grompp ion runfile

            # genion

            # translate box in z-dir

            # convert peptide coord file from gro to pdb

            # edit files prior to combining
            pass
            #

        # todo: implement algorithm here or actually in algorithm class?
        # todo: consider what is stored as the final element in a history list - indicate what is currently happening with the thread

        return running # todo: what will jobtype.algorithm return?

    def gatekeeper(self, thread, allthreads, settings):
        # todo: consider adding a function w/in gatekeeper that can identify an adsorption event and bring about early termination
        # todo: may not be necessary if it is decided that nvt will only be running for a short amount of time anyway

        # If job for this thread has status 'C'ompleted/'C'anceled...
        if thread.get_status(0, settings) == 'C':  # index 0 because there is only ever one element in thread.jobids
            # todo: implement restarting if simulation crashed before a certain number of steps completed?
            return True
        else:
            return False

