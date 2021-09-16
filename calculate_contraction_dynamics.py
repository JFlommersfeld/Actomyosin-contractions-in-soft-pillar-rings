import sys
from scipy.integrate import solve_ivp
from models import FullModel, DensityModel
from work_and_power import *


def calculate_contraction_dynamics(model_type, parameter_file, t_max, pillar_stiffness):
    """ Solves the contraction dynamics from time t=0 to t=t_max.

    Parameters:
        model_type (str):
            sets the type of model, which determines the exact parameter set that is needed. Possible values for
            the parameter model_type are: 'full model' or 'density model'.
        parameter_file (str):
            the path to the parameter file.
        t_max (float):
            the final timepoint in seconds until which the dynamics are solved.
        pillar_stiffness (float)
            the stiffness of the pillars in the pillar ring in pN/um.
    Returns:
        numpy.array(float): the calculated time points
        numpy.array(float): the solution of the (set of) differential equation(s)
        list(float): the tip velocities
        list(float): the transmitted powers
        list(float): the dissipated powers
        float: the total transmitted work
        float: the total dissipated work
    """
    if model_type == 'full model':
        model = FullModel(parameter_file, pillar_stiffness)
    elif model_type == 'density model':
        model = DensityModel(parameter_file, pillar_stiffness)
    else:
        print("ERROR: the parameter <model_type> has to be given one of the three values: "
              "'full model' or 'density model'.")
        sys.exit(1)

    t_min = 0

    # define initial conditions
    f_init = 7
    if model_type == 'full model':
        k_on_fil = model.get_parameter('k_on_fil')
        Nmax = model.get_parameter('Nmax')
        ic = [f_init, k_on_fil / (model.k_off_fil(f_init) + k_on_fil) * Nmax]
    else:
        ic = [0.]

    # integrate ODE
    sol = solve_ivp(model.rhs, (t_min, t_max), ic, method="LSODA")

    # calculate the dynamic quantities
    if model_type == 'full model':
        velocities = model.velocity(sol.t, sol.y[0], sol.y[1])
    else:
        velocities = model.velocity(sol.t, sol.y[0])

    transmitted_powers = transmitted_power(sol.y[0]/pillar_stiffness, velocities, model)
    dissipated_powers = dissipated_power(velocities, model)

    # calculate the total transmitted and dissipated work
    transmitted_work = work_elastic(sol.t, sol.y[0]/pillar_stiffness, velocities, model)
    dissipated_work = work_diss(sol.t, velocities, model)

    return sol.t, sol.y, velocities, transmitted_powers, dissipated_powers, transmitted_work, dissipated_work
