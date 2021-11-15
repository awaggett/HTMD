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
import interpret, process, jobtype
from infrastructure import configure, factory

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
        self.name = ''                  # name of current step
        self.jobids = []                # list of jobids associated with the present step of this thread
        self.terminated = False         # boolean indicating whether the thread has reached a termination criterion
        self.init_coords = ''
        self.peptides = []              # list of sequences (list of str) assigned to this thread for processing
        self.current_type = ''          # str indicating job type for the present step of this thread ('petide' or 'system')
        self.current_peptide = 0                # integer indicating index of current peptide in self.peptides

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
        return jobtype.get_batch_template(self, settings)

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
    peptide_groups = [pep_array.tolist() for pep_array in [x for x in np.array_split(settings.peptides, settings.nodes) if x.size > 0]]

    # number of peptide groups is equal to number of nodes
    for i, peptide_group in enumerate(peptide_groups):     # todo: this is arbitrary, for illustration purposes
        thread = Thread()
        #thread.init_coords = ''     # todo: assign initial coordinates (and other desired parameters) to each thread
        thread.peptides = peptide_group
        thread.name = settings.name + '_' + str(i)
        thread.current_peptide = 0                                  # index of current peptide in thread.peptides
        thread.current_type = ''  # todo: can maybe just start with initializing the current type as 'peptide'?
        thread.current_struct = 0  # not used by Adsorption_Simple class
        thread.current_rep = 0   # not used by Adsorption_Simple class
        jobtype.update_history(thread, settings, **{'initialize': True, 'add_peptides': peptide_group}) # todo: determine if any other kwargs needed

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

    # todo: moving a bunch of files for now, but find out a way to make this not necessary     
    # Move ff to working directory 
    if not os.path.exists(settings.working_directory + '/' + settings.force_field):
        source = settings.path_to_input_files + '/' + settings.force_field
        destination = settings.working_directory + '/' + settings.force_field 
        shutil.copytree(source, destination)

    if not os.path.exists(settings.working_directory + '/ion.mdp'):
        source = settings.path_to_input_files + '/ion.mdp'
        destination = settings.working_directory + '/ion.mdp'
        shutil.copyfile(source, destination)

    if not os.path.exists(settings.working_directory + '/QF_S_Silica_ph7.5_part.gro'):
       source = settings.path_to_input_files + '/QF_S_Silica_ph7.5_part.gro'
       destination = settings.working_directory + '/QF_S_Silica_ph7.5_part.gro'
       shutil.copyfile(source, destination)
     

    # Build or load threads
    allthreads = init_threads(settings)
    
    # Move runtime to working directory
    os.chdir(settings.working_directory) # todo: might I want different directories for peptides?
    # todo: should there be a check to make sure the forcefield is in the working directory?

    running = allthreads.copy()     # to be pruned later by thread.process()
    termination_criterion = False   # initialize global termination criterion boolean
    jobtype = factory.jobtype_factory(settings.job_type)    # initialize jobtype
    
    # Initialize threads with first process step
    try:
        print('FIRST STEP')
        for thread in allthreads:
            # todo: change if not thread.history.traj (instead check if thread.history.peptides matches original thread.peptides?)
            #if not thread.history.trajs:    # if there have been no steps in this thread yet
            if not any(thread.history.coords):   # if there have been no steps in this thread yet
                print('initializing history for: ', thread.peptides)
                jobtype.algorithm(thread, allthreads, running, settings)
            print('Ready to process: ', thread.peptides)
            running = thread.process(running, allthreads, settings)
    except Exception as e:
        if settings.restart:
            print('The following error occurred while attempting to initialize threads from restart.pkl. It may be '
                  'corrupted.')
        raise e

    return None
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
        print('working dir is: ' + str(sys.argv[2]))
    except IndexError:
        working_directory = ''
    print('config file is: ' + str(sys.argv[1]))
    settings = configure.configure(sys.argv[1], working_directory)
    exit_message = main(settings)
    print(exit_message)

if __name__ == "__main__":
    run_main()

