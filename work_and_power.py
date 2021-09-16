import numpy as np
from math import pi


def transmitted_power(displacements, velocities, model):
    """Calculates the power that is transmitted to the pillars

    Parameters:
        displacements (numpy.array(float)):
            list of the displacements of the pillar tip.
        velocities (numpy.array(float)):
            list of the velocities of the pillar tip.
        model:
            a model object that implements the get_parameter function.

    Returns:
        float: the transmitted power
    """
    k_p = model.get_parameter('k_p')
    powers = np.asarray(5. * 2. * k_p * displacements * velocities)
    return powers


def dissipated_power(velocities, model):
    """Calculates the power that is dissipated by the viscous gel

    Parameters:
        velocities (numpy.array(float)):
            list of the velocities of the pillar tip.
        model:
            a model object that implements the get_parameter function.
    Returns:
        float: the transmitted power
    """
    gamma = 6. / 5. * pi * model.get_parameter('h_eta')
    powers = np.asarray(5. * gamma * 4. * velocities ** 2)
    return powers


def work_elastic(t, displacements, velocities, model):
    """Calculates the work that is performed on the pillars throughout a given time series

    Parameters:
        t (numpy.array(float)):
            list of timepoints
        displacements (numpy.array(float)):
            list of the displacements of the pillar tip.
        velocities (numpy.array(float)):
            list of the velocities of the pillar tip.
        model:
            a model object that implements the get_parameter function.

    Returns:
        float: the work that was perfomred on the pillars
    """
    work = 0
    powers = transmitted_power(displacements, velocities, model)
    time_prev = 0
    for (time, power) in zip(t, powers):
        work += (time - time_prev) * power
    return work


def work_diss(t, velocities, model):
    """Calculates the work that is dissipated by the viscous gel throughout a given time series

    Parameters:
        t (numpy.array(float)):
            list of timepoints
        velocities (numpy.array(float)):
            list of the velocities of the pillar tip.
        model:
            a model object that implements the get_parameter function.
    Returns:
        float: the work that was dissipated by the viscous gel
    """
    work = 0
    powers = dissipated_power(velocities, model)
    time_prev = 0
    for (time, power) in zip(t, powers):
        work += (time - time_prev) * power
    return work
