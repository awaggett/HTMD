"""
Interface for Algorithm objects. New Algorithms can be implemented by constructing a new class that inherits from
Algorithm and implements its abstract methods.
"""

import os
import re
import sys
import abc
import copy
import time
import numpy
import shutil
import pickle
import mdtraj
import argparse
import itertools
import subprocess
import utilities
from math import factorial
#from filelock import Timeout, FileLock

class Algorithm(abc.ABC):
    """
    Abstract base class for isEE algorithms.
    Implements methods for all of the algorithm-specific tasks that isEE might need.
    """

    def __init__(self):
        pass

    @staticmethod
    def build_algorithm_history(allthreads):
        # First, establish lock so only one instance of build_algorithm_history will run at a time
        lock = FileLock('algorithm_history.lock')
        with lock:
            open('arbitrary_lockfile.lock', 'w').close()
        lock.acquire()
        open('arbitrary_lockfile.lock', 'w').close()    # this step will block until the lock is released

        # Build algorithm_history afresh from all thread history attributes
        if settings.shared_history_file:
            history_file = settings.shared_history_file
        else:
            history_file = 'algorithm_history.pkl'

        if not os.path.exists(history_file):
            algorithm_history = argparse.Namespace()
        else:
            algorithm_history = pickle.load(open(history_file, 'rb'))
        # todo: check that this still applies


        def merge_namespace(current, new):
            # Helper function to merge contents of 'new' namespace into 'current' namespace, assuming both have the same
            # keys, that each key corresponds to a list, and that there are the same number of entries in each list.
            # If a particular "column" of entries in 'new' is already present in 'current' then it is not merged.]))
            # 'current' is modified in place so there is no need for a return statement
            for key in list(current.__dict__.keys()):
                current.__dict__[key] += [new.__dict__[key][col]                                                # append dictionary entry with key 'key' in column 'col' ...
                                          for col in range(len(new.__dict__[list(new.__dict__.keys())[0]]))     # for all columns 'col' in 'new' ...
                                          if not all([new.__dict__[kk][col] in current.__dict__[kk]             # if this column does not already match the contents of any other column ...
                                                      for kk in list(current.__dict__.keys())])]                # for all keys

        # Perform merge
        for thread in allthreads:
            merge_namespace(algorithm_history, thread.history)

        # Sort every column chronologically by timestamp attribute
        for key in list(algorithm_history.__dict__.keys()):
            algorithm_history.__dict__[key] = [x for _,x in sorted(zip(algorithm_history.__dict__['timestamps'], algorithm_history.__dict__[key]))]

        # Dump pickle file to 'algorithm_history.pkl' regardless of settings.shared_history_file, because every working
        # directory should have its own algorithm history object in addition to the shared one.
        pickle.dump(algorithm_history, open('algorithm_history.pkl.bak', 'wb'))
        if not os.path.getsize('algorithm_history.pkl.bak') == 0:
            shutil.copy('algorithm_history.pkl.bak', 'algorithm_history.pkl')  # copying after is safer
        else:
            raise RuntimeError('failed to dump algorithm history pickle file')

        # Copy to shared file if necessary
        if settings.shared_history_file:
            shutil.copy('algorithm_history.pkl', settings.shared_history_file)

        lock.release()
        return algorithm_history
        # todo: check that this still applies

    @abc.abstractmethod
    def get_first_step(self, thread, allthreads, settings):
        """
        Determine the first step for a thread. This method should be called before the first 'process' step after a new
        thread is created (as defined by having an empty thread.history.trajs attribute) to allow the algorithm to
        decide what its first step should be.
        Specifically, implementations of this method should either pass on all threads after the first one directly to
        get_next_step (for algorithms where the next step can be determined without information from the first one), or
        else idle each thread after the first one.
        Parameters
        ----------
        thread : Thread()
            Thread object to consider
        allthreads : list
            List of all thread object to consider
        settings : argparse.Namespace
            Settings namespace object
        Returns
        -------
        first_step : list or str

        """

        pass

    @abc.abstractmethod
    def get_next_step(self, thread, allthreads, settings):
        """
        Determine the next step for the thread,

        Parameters
        ----------
        thread : Thread()
            Thread object to consider
        allthreads : list
            Full list of Thread() objects for this job (including thread)
        settings : argparse.Namespace
            Settings namespace object
        Returns
        -------
        next_step : list or str

        """

        pass


class Adsorption(Algorithm):
    """
    Adapter class for algorithm that:
        - Builds peptide using Amber TEaP
        - Minimizes and equilibrates peptide using Gromacs
        - Grafts relaxed peptide onto user's choice of surface
        - Minimizes and equilibrates peptide-surface system
    """

    def get_first_step(self, thread, settings):
        # build peptide according to sequence. utilities.build_peptide
        # edit amber pdb for gromacs. utilities.edit_pdb
        pass
