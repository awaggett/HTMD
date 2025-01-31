"""
This portion of the program is responsible for handling update of the results, checking global termination criteria, and
implementing the calls to JobType methods to control the values of the thread attributes for the next step.
"""

import os
import shutil
import pickle
from infrastructure import factory

def interpret(thread, allthreads, running, settings):
    """
    The main function of interpret.py. Makes calls to JobType methods to update results, check termination criteria, and
    update thread parameters

    Parameters
    ----------
    thread : Thread
        The Thread object on which to act
    allthreads : list
        The list of all extant Thread objects
    running : list
        The list of all currently running Thread objects
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    termination: bool
        True if a global termination criterion has been met; False otherwise

    """

    jobtype = factory.jobtype_factory(settings.job_type)

    jobtype.analyze(thread, settings)                               # analyze just-completed job step
    termination = jobtype.algorithm(thread, allthreads, running, settings)   # query algorithm to decide next move (or terminate)
    # todo: in the above step, it is going to move the current step to the next and build the next system (will need to shift thread.current_peptide and other props there)
    return termination, running
