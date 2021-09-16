import matplotlib.pyplot as plt


def plot_tip_displacement_and_velocity(timepoints, displacements, velocities, ylim_bottom_velo=0, ylim_top_velo=10,
                                       ylim_bottom_disp=0, ylim_top_disp=7):
    """
    Plots the dynammics of the tip displacement and tip velocity

        Parameters:
            timepoints (list(float)):
                list of the considered timepoints in seconds.
            displacements (list(float)):
                list of the calculated tip displacements in um.
            velocities (list(float)):
                list of the calculated tip velocities in um/s.
            ylim_bottom_velo (float):
                the lower limit of the y-axis for the tip velocity axis.
            ylim_top_velo (float):
                the upper limit of the y-axis for the tip velocity axis.
            ylim_bottom_disp (float):
                the lower limit of the y-axis for the tip displacment axis.
            ylim_top_disp (float):
                the upper limit of the y-axis for the tip displacment axis.
    """
    plt.rcParams.update({'font.size': 15})

    plt.clf()
    fig, ax1 = plt.subplots()

    color1 = 'k'
    ax1.set_xlabel('Time (min)')
    ax1.set_ylabel(r'Tip velocity ($10^{-2}\mu m /s$)', color=color1)
    ax1.plot(timepoints/60, 100*velocities, color=color1, linewidth=2)

    ax1.tick_params(axis='y', labelcolor=color1, color=color1)
    ax1.set_ylim(bottom=ylim_bottom_velo, top=ylim_top_velo)

    ax2 = ax1.twinx()

    color2 = 'tab:red'
    ax2.set_ylabel(r'Tip displacement ($\mu m$)', color=color2)
    ax2.plot(timepoints/60, displacements, color=color2, linewidth=2)
    ax2.tick_params(axis='y', labelcolor=color2, color=color2)
    ax2.spines['left'].set_color(color1)
    ax2.spines['right'].set_color(color2)
    ax2.set_ylim(bottom=ylim_bottom_disp, top=ylim_top_disp)

    fig.tight_layout()
    plt.show()


def plot_transmitted_and_dissipated_power(timepoints, transmitted_powers, dissipated_powers, ylim_bottom=0.1,
                                          ylim_top=1.e3):
    """
    Plots the dynammics of the transmitted and dissipated power on a logarithmic y-axis.

        Parameters:
            timepoints (list(float)):
                list of the considered timepoints in seconds.
            transmitted_powers (list(float)):
                list of the calculated transmitted powers in atto Watts.
            dissipated_powers (list(float)):
                list of the calculated dissipated powers atto Watts.
            ylim_bottom (float):
                the lower limit of the y-axis. Default value is 0.1.
            ylim_top (float):
                the upper limit of the y-axis. Default value is 1.e3.
    """
    plt.rcParams.update({'font.size': 15})

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Time (min)')

    ax1.set_ylabel(r'Power ($aW$)')
    ax1.semilogy(timepoints/60, transmitted_powers, 'k-', linewidth=2, alpha=1, label='Transmitted power')
    ax1.semilogy(timepoints/60, dissipated_powers, 'k:', linewidth=2, alpha=1, label='Dissipated power')
    ax1.legend(frameon=False, loc=(0.5, 0.75), fontsize=13)
    ax1.set_ylim(bottom=ylim_bottom, top=ylim_top)

    fig.tight_layout()
    plt.show()


def plot_peak_velocities(pillar_stiffnesses, peak_velocities, ylim_bottom=0., ylim_top=10.):
    """
    Plots the dynammics of the transmitted and dissipated power

        Parameters:
            pillar_stiffnesses (list(float)):
                list of the considered pillar stiffnesses in unites of pN/um.
            peak_velocities (list(float)):
                list of the calculated peak velocities in um/s.
            ylim_bottom (float):
                the lower limit of the y-axis. Default value is 0.
            ylim_top (float):
                the upper limit of the y-axis. Default value is 10.
    """
    plt.plot(pillar_stiffnesses, 100 * peak_velocities, '-', linewidth=2)
    plt.xlabel(r"Pillar stiffness $k_p$ ($pN/\mu m$)")
    plt.ylabel(r"Peak velocity  ($10^{-2}\mu m/s $)")
    plt.ylim(bottom=ylim_bottom, top=ylim_top)
    plt.tight_layout()
    plt.show()


def plot_final_forces(pillar_stiffnesses, final_force, ylim_bottom=0., ylim_top=500.):
    """
    Plots the dynammics of the transmitted and dissipated power

        Parameters:
            pillar_stiffnesses (list(float)):
                list of the considered pillar stiffnesses in unites of pN/um.
            final_force (list(float)):
                list of the calculated final forces per pillar in pN.
            ylim_bottom (float):
                the lower limit of the y-axis. Default value is 0.
            ylim_top (float):
                the upper limit of the y-axis. Default value is 500.
    """
    plt.clf()
    plt.rcParams.update({'font.size': 15})
    plt.plot(pillar_stiffnesses, final_force)
    plt.xticks([50, 100, 150], ["0.05", "0.10", "0.15"])

    plt.xlim(left=0, right=190)
    plt.ylim(ylim_bottom, ylim_top)
    plt.xlabel(r"Pillar stiffness $k_p$ ($nN/\mu m$)")
    plt.ylabel(r"Force per pillar ($pN$)")

    plt.tight_layout()
    plt.show()
