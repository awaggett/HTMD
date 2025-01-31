"""
Factory script for obtaining the desired interfaces from the various interface scripts.
"""

import jobtype
import algorithm
#from htmd import stabilitymodel
from infrastructure import mdengine
from infrastructure import batchsystem
from infrastructure import taskmanager

def mdengine_factory(mdengine_toolkit):
    """
    Factory function for MDEngines.
    Parameters
    ----------
    mdengine_toolkit : str
        Name of the MDEngine to invoke
    Returns
    -------
    mdengine : MDEngine
        Instance of an MDEngine adapter
    """

    mdengine_toolkit = mdengine_toolkit.lower()

    mdengine_toolkits = {'gromacs': mdengine.AdaptGromacs()} # todo: do I need an MDengine?

    if mdengine_toolkit not in mdengine_toolkits.keys():
        raise ValueError('unsupported MDEngine name: ' + mdengine_toolkit)

    return mdengine_toolkits[mdengine_toolkit]


def batchsystem_factory(batchsystem_toolkit):
    """
    Factory function for BatchSystems.
    Parameters
    ----------
    batchsystem_toolkit : str
        Name of the BatchSystem to invoke
    Returns
    -------
    batchsystem : BatchSystem
        Instance of a BatchSystem adapter
    """

    batchsystem_toolkit = batchsystem_toolkit.lower()

    batchsystem_toolkits = {'slurm': batchsystem.AdaptSlurm(),
                            'pbs': batchsystem.AdaptPBS(),
                            'torque': batchsystem.AdaptPBS()}   # torque and pbs are synonyms

    if batchsystem_toolkit not in batchsystem_toolkits.keys():
        raise ValueError('unsupported BatchSystem name: ' + batchsystem_toolkit)

    return batchsystem_toolkits[batchsystem_toolkit]


def jobtype_factory(jobtype_toolkit):
    """
    Factory function for JobTypes.
    Parameters
    ----------
    jobtype_toolkit : str
        Name of the JobType to invoke
    Returns
    -------
    jobtype : JobType
        Instance of a JobType adapter
    """

    jobtype_toolkit = jobtype_toolkit.lower()

    jobtype_toolkits = {'adsorption': jobtype.Adsorption()}

    if jobtype_toolkit not in jobtype_toolkits.keys():
        raise ValueError('unsupported JobType name: ' + jobtype_toolkit)

    return jobtype_toolkits[jobtype_toolkit]

def taskmanager_factory(taskmanager_toolkit):
    """
    Factory function for TaskManagers.
    Parameters
    ----------
    taskmanager_toolkit : str
        Name of the TaskManager to invoke
    Returns
    -------
    taskmanager : TaskManager
        Instance of a TaskManager adapter
    """

    taskmanager_toolkit = taskmanager_toolkit.lower()

    taskmanager_toolkits = {'simple': taskmanager.AdaptSimple()}

    if taskmanager_toolkit not in taskmanager_toolkits.keys():
        raise ValueError('unsupported TaskManager name: ' + taskmanager_toolkit)

    return taskmanager_toolkits[taskmanager_toolkit]

def algorithm_factory(algorithm_toolkit):
    """
    Factory function for Algorithms.
    Parameters
    ----------
    algorithm_toolkit : str
        Name of the Algorithm to invoke
    Returns
    -------
    algorithm : Algorithm
        Instance of a Algorithm adapter
    """

    algorithm_toolkit = algorithm_toolkit.lower()

    algorithm_toolkits = {'script': algorithm.Script(),
                          'covariance_saturation': algorithm.CovarianceSaturation(),
                          'subnetwork_hotspots': algorithm.SubnetworkHotspots(),
                          'monte_carlo': algorithm.MonteCarlo()}

    if algorithm_toolkit not in algorithm_toolkits.keys():
        raise ValueError('unsupported Algorithm name: ' + algorithm_toolkit)

    return algorithm_toolkits[algorithm_toolkit]


