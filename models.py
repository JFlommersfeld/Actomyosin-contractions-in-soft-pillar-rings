import numpy as np
from parameter_loader import load_parameters
from math import pi


class FullModel:
    """
    A class that defines the full model for actomyosin contraction in soft pillar rings that accounts for both myosin
    filament binding and density changes

    Attributes:
        parameter_file (str):
            path to a file that contains all necessary parameters for the model (see provided examples).

    Methods:
        k_off_fil(total_force):
            calculates the load dependent steady state off-rate of a myosin filament.

        rhs(t, y):
            calculates the right hand side of the set of differential equations that describe the model.
    
        velocity(t, force, N):
            calculates the deflection velocity of the tip of the pillar.
    """
    
    def __init__(self, parameter_file, pillar_stiffness):
        """
        Sets all the necessary parameters for the FullModel object.

        Parameters:
            parameter_file (str):
                path to a file that contains all necessary parameters for the model (see provided examples).
            pillar_stiffness (float):
                stiffness of the pillars in the pillar ring in pN/um.
        """
        self.x_catch, self.x_slip, self.k_off0_catch, self.k_off0_slip, self.k_on, self.k_on_fil, self.a_per_kBT, \
            self.Nh, self.Nmax, self.h_eta, self.xi_rho_a2, self.rho_max_per_rho, \
            self.R0 = load_parameters('full model', parameter_file)
        self.k_p = pillar_stiffness

        self.parameter_dict = {"x_catch": self.x_catch, "x_slip": self.x_slip, "k_off0_catch": self.k_off0_catch,
                               "k_off0_slip": self.k_off0_slip, "k_on": self.k_on, "k_on_fil": self.k_on_fil,
                               "a_per_kBT": self.a_per_kBT, "Nh": self.Nh, "Nmax": self.Nmax, "h_eta": self.h_eta,
                               "xi_rho_a2": self.xi_rho_a2, "rho_max_per_rho": self.rho_max_per_rho, "R0": self.R0,
                               "k_p": self.k_p}

        self.A0 = pi * self.R0**2
        self.tau = 6. / 5. * pi * self.h_eta / self.k_p

    def __k_off(self, force):
        """Calculates the load dependent off-rate of an individual myosin head.
        
        Parameters:
            force (float):
                the average load that is applied to an individual myosin head.

        Returns:
            float: the average off-rate of the head.
        """
        return self.k_off0_catch * np.exp(-self.a_per_kBT * force * self.x_catch) + \
            self.k_off0_slip * np.exp(self.a_per_kBT * force * self.x_slip)

    def __calc_prob_dist(self, total_force):
        """Calculates the load dependent steady state probability distribution of the number of bound heads per
           myosin filament
        
        Parameters:
            total_force (float):
                the total load that is applied to the myosin filament.

        Returns:
            list(float): list of probabilities that n heads are bound per filament, where n is given by the list index.
        """
        pns = []
        for n in range(0, self.Nh + 1):
            nom = 1
            for i in range(0, n):
                nom = nom * ((self.Nh - i) * self.k_on) / ((i + 1) * self.__k_off(total_force / (i + 1)))
            
            denom = 1
            for k in range(1, self.Nh + 1):
                prod = 1
                for j in range(0, k):
                    prod = prod * ((self.Nh - j) * self.k_on) / ((j + 1) * self.__k_off(total_force / (j + 1)))
                denom = denom + prod
            
            pns.append(nom / denom)
        
        return pns
        
    def k_off_fil(self, total_force):
        """Calculates the load dependent steady state off-rate of a myosin filament.
        
        Parameters:
            total_force (float):
                the total load that is applied to the myosin filament.

        Returns:
            float: the off-rate of the filament.
        """
        T_off_av = 0
        pns = self.__calc_prob_dist(total_force)
        for NB_init in range(1, self.Nh + 1):
            T_off = 0
            for NB in range(1, NB_init + 1):
                s = 0
                for j in range(NB, self.Nh + 1):
                    s = s + pns[j]
                
                T_off = T_off + 1 / (NB * self.__k_off(total_force / NB) * pns[NB]) * s
            
            T_off_av = T_off_av + pns[NB_init] * T_off
        return 1 / T_off_av

    def rhs(self, t, y):
        """Calculates the right hand side of the set of differential equations that describe the model.
        
        Parameters:
            t (float):
                the time point.
            y (list(float)):
                a list with elements y[0] = force on the pillar at time t and y[1] = number of bound filaments at time t

        Returns:
            list(float): the temporal derivative of the input y
        """
        force = y[0]
        N = y[1]

        area = pi * (self.R0 - force / self.k_p) ** 2
        density_factor = -self.A0 / area * (self.A0 / area - self.rho_max_per_rho)

        force_prime = -force / self.tau + self.xi_rho_a2 * N * density_factor / self.tau
        N_prime = self.k_on_fil * (self.Nmax - N) - self.k_off_fil(force) * N

        return [force_prime, N_prime]

    def velocity(self, t, force, N):
        """Calculates the deflection velocity of the tip of the pillar.
        
        Parameters:
            t (float):
                the time point.
            force (float):
                force on the pillar at time t
            N:
                number of bound filaments at time t

        Returns:
            float: the deflection velocity of the pillar tip at time t
        """
        area = pi * (self.R0 - force / self.k_p) ** 2
        density_factor = -self.A0 / area * (self.A0 / area - self.rho_max_per_rho)

        return (-force / self.tau + self.xi_rho_a2 * N * density_factor / self.tau) / self.k_p

    def get_parameter(self, parameter_name):
        """Get all model parameters

        Parameters:
            parameter_name (str):
                parameter name.

        Returns:
            float/int: the value of the specified parameter.
        """

        return self.parameter_dict[parameter_name]


class DensityModel:
    """
    A class that defines the purley density dependent model for actomyosin contraction in soft pillar rings.
    ...

    Attributes:
        parameter_file (str):
            path to a file that contains all necessary parameters for the model (see provided examples).

    Methods:
        k_off_fil(total_force):
            calculates the load dependent steady state off-rate of a myosin filament.

        rhs(t, y):
            calculates the right hand side of the set of differential equations that describe the model.

        velocity(t, force, N):
            calculates the deflection velocity of the tip of the pillar.
    """
    
    def __init__(self, parameter_file, pillar_stiffness):
        """
        Sets all the necessary parameters for the DensityModel object.

        Parameters:
            parameter_file (str):
                path to a file that contains all necessary parameters for the model (see provided examples).
            pillar_stiffness (float):
                stiffness of the pillars in the pillar ring in pN/um.
        """
        self.h_eta, self.xi_N_rho_a2, self.rho_max_per_rho, self.R0 = load_parameters('density model', parameter_file)
        self.k_p = pillar_stiffness

        self.parameter_dict = {"h_eta": self.h_eta, "xi_N_rho_a2": self.xi_N_rho_a2,
                               "rho_max_per_rho": self.rho_max_per_rho, "R0": self.R0, "k_p": self.k_p}

        self.A0 = pi * self.R0 ** 2
        self.tau = 6. / 5. * pi * self.h_eta / self.k_p

    def rhs(self, t, y):
        """Calculates the right hand side of the set of differential equations that describe the model.
        
        Parameters:
            t (float):
                the time point.
            y (list(float)):
                a list with a single element y[0] = force on the pillar at time t

        Returns:
        list(float): the temporal derivative of the input y
        """
        force = y[0]

        area = pi * (self.R0 - force / self.k_p) ** 2
        density_factor = -self.A0 / area * (self.A0 / area - self.rho_max_per_rho)
        
        force_prime = -force/self.tau + self.xi_N_rho_a2 * density_factor / self.tau

        return [force_prime]

    def velocity(self, t, force):
        """Calculates the deflection velocity of the tip of the pillar.
        
        Parameters:
            t (float):
                the time point.
            force (float):
                force on the pillar at time t

        Returns:
        float: the deflection velocity of the pillar tip at time t
        """
        area = pi * (self.R0 - force / self.k_p) ** 2
        density_factor = -self.A0 / area * (self.A0 / area - self.rho_max_per_rho)

        return (-force/self.tau + self.xi_N_rho_a2 * density_factor / self.tau)/self.k_p

    def get_parameter(self, parameter_name):
        """Get all model parameters

        Parameters:
            parameter_name (str):
                parameter name.

        Returns:
            float/int: the value of the specified parameter.
        """

        return self.parameter_dict[parameter_name]
