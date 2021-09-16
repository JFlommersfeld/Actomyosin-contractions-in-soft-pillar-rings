import numpy as np
from calculate_contraction_dynamics import calculate_contraction_dynamics
from plotting import plot_tip_displacement_and_velocity, plot_peak_velocities, plot_final_forces
from auxiliaries import discretize_curve, velos_numerical


# first we analyse the contraction dynamics for a single pillar stiffness
pillar_stiffness = 35.

# simulate the model
timepoints, sol, velocities, _, _, _, _ = \
    calculate_contraction_dynamics('density model', './parameter_density_model.txt', 210., pillar_stiffness)

# plot the dynamic quantities
plot_tip_displacement_and_velocity(timepoints, sol[0]/pillar_stiffness, velocities, ylim_bottom_velo=-0.5,
                                   ylim_top_velo=5, ylim_bottom_disp=-0.5, ylim_top_disp=4)

# repeat the analysis at a higher pillar stiffness
pillar_stiffness = 100.

# simulate the model
timepoints, sol, velocities, _, _, _, _ = \
    calculate_contraction_dynamics('density model', './parameter_density_model.txt', 210., pillar_stiffness)

# plot the dynamic quantities
plot_tip_displacement_and_velocity(timepoints, sol[0]/pillar_stiffness, velocities, ylim_bottom_velo=-0.5,
                                   ylim_top_velo=5, ylim_bottom_disp=-0.5, ylim_top_disp=4)


# now we analyse the dependence on the pillar stiffness
peak_velo_discrete = []
final_forces = []
t_max = 500.
for idx, stiffness in enumerate(range(20, 230, 5)):
    pillar_stiffness = stiffness

    # integrate ODE
    timepoints, sol, velocities, transmitted_powers, dissipated_powers, transmitted_work, dissipated_work = \
        calculate_contraction_dynamics('density model', './parameter_density_model.txt', t_max, pillar_stiffness)

    # get the final force and discretized peak velocity
    final_forces.append([pillar_stiffness, sol[0][-1]])
    times_discrete, disp_discrete = discretize_curve(timepoints, sol[0] / pillar_stiffness, 30)
    velos_discrete = velos_numerical(np.asarray(times_discrete) / 60, disp_discrete)
    peak_velo_discrete.append([pillar_stiffness, np.max(velos_discrete)])

# plot peak velocities
peak_velo_discrete = np.transpose(np.asarray(peak_velo_discrete))
plot_peak_velocities(peak_velo_discrete[0], peak_velo_discrete[1])

# plot final forces
final_forces = np.asarray(final_forces)
final_forces = final_forces.T
plot_final_forces(final_forces[0], final_forces[1], ylim_bottom=-20, ylim_top=420)
