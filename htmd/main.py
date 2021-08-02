"""
main.py
Non-functional framework for managing jobs in isEE/atesa

In comments and job structure throughout 'isEE' or 'htmd' will appear; this is a hold-over from the original isEE code.
THIS IS NOT isEE!

This script handles the primary loop of building and submitting jobs in independent Threads, using the methods thereof 
to execute various interfaced/abstracted commands.
"""

import copy
import os
import pickle
import shutil
import sys
import time
import numpy as np
from htmd import interpret, process, jobtype
from htmd.infrastructure import configure, factory

class Thread(object):
    """
    Object representing a series of simulations and containing the relevant information to define its current state.

    Threads represent the level on which isEE is parallelized. This flexible object is used for every type of job
    performed by isEE.

    Parameters
    ----------
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    None

    """

    def __init__(self):
        self.name = ''
        self.jobids = []                # list of jobids associated with the present step of this thread
        self.terminated = False         # boolean indicating whether the thread has reached a termination criterion
        self.init_coords = ''
        self.sequences = []             # list of sequences (list of str) assigned to this thread for processing
        self.current_sequence = ''      # todo: may change this but woudl be good for thread to track what it is currently processing

    # Remember in implementations of Thread methods that 'self' is a thread object, even though they may make calls to
    # methods of other types of objects (that therefore use 'self' to refer to non-Thread objects)

    def process(self, running, allthreads, settings):
        return process.process(self, running, allthreads, settings)

    def interpret(self, allthreads, running, settings):
        return interpret.interpret(self, allthreads, running, settings)

    def gatekeeper(self, allthreads, settings):
        jobtype = factory.jobtype_factory(settings.job_type)
        return jobtype.gatekeeper(self, allthreads, settings)

    def get_next_step(self, settings):
        jobtype = factory.jobtype_factory(settings.job_type)
        self.current_name = jobtype.get_next_step(self, settings)
        return self.current_name

    def get_batch_template(self, settings):
        jobtype = factory.jobtype_factory(settings.job_type)
        return jobtype.get_batch_template(settings)

    def get_frame(self, traj, frame, settings):
        mdengine = factory.mdengine_factory(settings.md_engine)
        return mdengine.get_frame(traj, frame, settings)

    def get_status(self, job_index, settings):
        batchsystem = factory.batchsystem_factory(settings.batch_system)
        return batchsystem.get_status(self.jobids[job_index], settings)

    def cancel_job(self, job_index, settings):
        batchsystem = factory.batchsystem_factory(settings.batch_system)
        batchsystem.cancel_job(self.jobids[job_index], settings)


def init_threads(settings):
    """
    Initialize all the Thread objects called for by the user input file.

    In the case where settings.restart == True, this involves unpickling restart.pkl; otherwise, brand new objects are
    produced in accordance with settings.job_type (aimless_shooting, committor_analysis, equilibrium_path_sampling, or
    htmd).

    Parameters
    ----------
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    allthreads : list
        List of all Thread objects, including those for which no further tasks are scheduled.

    """

    if settings.restart:
        allthreads = pickle.load(open(settings.working_directory + '/restart.pkl', 'rb'))

        # todo: I like to save thread objects as pickle files every so often so that they can be reloaded like this
        # todo: in order to restart crashed or otherwise terminated runs. This may require some additional 
        # todo: processing on your part.
 
        return allthreads   # this return statement precludes the remainder of this function from running

    # If not restart:
    allthreads = []
    jobtype = factory.jobtype_factory(settings.job_type)

    # todo: here, you should initialize threads as desired using something like this:
    # partition peptides equally as possible between threads based on the number of nodes and peptides to be processed
    peptide_groups = [pep_array.tolist() for pep_array in np.array_split(settings.peptides,settings.nodes)]

    for i, peptide_group in enumerate(peptide_groups):     # todo: this is arbitrary, for illustration purposes
        thread = Thread()
        #thread.init_coords = ''     # todo: assign initial coordinates (and other desired parameters) to each thread
        thread.sequences = peptide_group
        thread.name = settings.name + '_' + str(i)
        allthreads.append(thread)

    return allthreads

# This function is used to cancel submitted jobs in the event that the parent python program crashes
def handle_loop_exception(running, exception, settings):
    """
    Handle cancellation of jobs after encountering an exception.

    Parameters
    ----------
    running : list
        List of Thread objects that are currently running. These are the threads that will be canceled if the isEE run
        cannot be rescued.
    exception : Exception
        The exception that triggered calling this function
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    None

    """

    print('\nCancelling currently running batch jobs belonging to this process in order to '
          'preserve resources.')
    for thread in running:
        try:
            for job_index in range(len(thread.jobids)):
                thread.cancel_job(job_index, settings)
        except Exception as little_e:
            print('\nEncountered an additional exception while attempting to cancel a job: ' + str(little_e) +
                  '\nIgnoring and continuing...')

    print('Job cancellation complete, isEE is now shutting down. The full exception that triggered this was: ')

    raise exception


def main(settings):
    """
    Perform the primary loop of building, submitting, monitoring, and analyzing jobs.

    This function works via a loop of calls to thread.process and thread.interpret for each thread that hasn't
    terminated, until either the global termination criterion is met or all the individual threads have completed.

    Parameters
    ----------
    settings : argparse.Namespace
        Settings namespace object
    rescue_running : list
        List of threads passed in from handle_loop_exception, containing running threads. If given, setup is skipped and
        the function proceeds directly to the main loop.

    Returns
    -------
    None

    """

    # Make working directory if it does not exist, handling overwrite and restart as needed
    if os.path.exists(settings.working_directory):
        if settings.overwrite and not settings.restart:
            shutil.rmtree(settings.working_directory)
            os.mkdir(settings.working_directory)
        elif not settings.restart:
            raise RuntimeError('Working directory ' + settings.working_directory + ' already exists, but overwrite '
                               '= False and restart = False. Either change one of these two settings or choose a '
                               'different working directory.')
    else:
        if not settings.restart:
            os.mkdir(settings.working_directory)
        else:
            raise RuntimeError('Working directory ' + settings.working_directory + ' does not yet exist, but '
                               'restart = True.')

    # Build or load threads
    allthreads = init_threads(settings)

    # Move runtime to working directory
    os.chdir(settings.working_directory)

    running = allthreads.copy()     # to be pruned later by thread.process()
    termination_criterion = False   # initialize global termination criterion boolean
    #jobtype = factory.jobtype_factory(settings.job_type)    # initialize jobtype
    jobtype = jobtype.Adsorption()

    # Initialize threads with first process step
    try:
        for thread in allthreads:
            if not thread.history.trajs:    # if there have been no steps in this thread yet
                jobtype.algorithm(thread, allthreads, settings)
            running = thread.process(running, allthreads, settings)
    except Exception as e:
        if settings.restart:
            print('The following error occurred while attempting to initialize threads from restart.pkl. It may be '
                  'corrupted.')
        raise e

    # Begin main loop
    # This whole thing is in a try-except block to handle cancellation of jobs when the code crashes in any way
    try:
        while (not termination_criterion) and running:
            for thread in running:
                if thread.gatekeeper(allthreads, settings):
                    termination_criterion, running = thread.interpret(allthreads, running, settings)
                    if termination_criterion:   # global termination
                        for thread in running:
                            for job_index in range(len(thread.jobids)):
                                thread.cancel_job(job_index, settings)
                        running = []
                        break
                    running = thread.process(running, allthreads, settings)
                else:
                    time.sleep(30)  # to prevent too-frequent calls to batch system by thread.gatekeeper

    except Exception as e:  # catch arbitrary exceptions in the main loop so we can cancel jobs
        print(str(e))
        handle_loop_exception(running, e, settings)

    if termination_criterion:
        return 'isEE run exiting normally (global termination criterion met)'
    else:
        return 'isEE run exiting normally (all threads ended individually)'


def run_main():
    # Obtain settings namespace, initialize threads, and move promptly into main.
    try:
        working_directory = sys.argv[2]
    except IndexError:
        working_directory = ''
    settings = configure.configure(sys.argv[1], working_directory)
    exit_message = main(settings)
    print(exit_message)

if __name__ == "__main__":
    run_main()
