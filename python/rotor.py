
import numpy as np
from scipy import optimize

class Rotor:
    """
    Rotors as aircraft components.

    Attributes
    ----------
    radius : float
        Rotor radius; m
    number_of_blades : int
        Number of blades.
    chord : float
        Chord length; m
    kappa : float
        Factor between ideal and real induced power.
    zero_lift_drag_coeff : float
        Zero-lift drag coefficient of the rotor blade, average.
    tip_velocity : float
        Tip velocity; m/s
    power_fraction : float
        Power fraction relative to the main rotor, if applicable.
    installation_height : float
        Distance between landing gear and rotor; m

    """
    
    def __init__(self, rotor_data: dict):
        """
        Load the rotor variables.

        Parameters
        ----------
        rotor_data : dict
            Rotor section of the configuration file.

        """
        self.radius               = rotor_data.get('Radius')
        self.number_of_blades     = rotor_data.get('Number of blades')
        self.chord                = rotor_data.get('Chord')
        self.kappa                = rotor_data.get('Kappa', 1.0)
        self.zero_lift_drag_coeff = rotor_data.get('Zero-lift drag coeff.')
        self.tip_velocity         = rotor_data.get('Tip velocity')
        self.power_fraction       = rotor_data.get('Power fraction')
        self.installation_height  = rotor_data.get('Installation height')


    @property
    def solidity(self):
        """
        Rotor solidity (rectangular approximation).

        Note
        ----
        The solidity is implemented as a property, which means it can be used
        as an attribute, calculated on call. This way it will always be in
        sync with the rotor radius and chord length.
        
        Returns
        -------
        float
            Solidity.

        """
        # Solidity [-]
        return self.number_of_blades * self.chord / (np.pi * self.radius)


    def get_disc_loading(self, thrust):
        """
        Calculate the disc loading.
        
        Parameters
        ----------
        thrust : float
            Thrust; N
        
        Returns
        -------
        float
            Disc loading; N/m²

        """
        # Disc loading [N/m²]
        return thrust / (np.pi * self.radius ** 2)


    def get_blade_loading(self, density, thrust):
        """
        Calculate the blade loading CT / sigma. [p.133]
        
        Parameters
        ----------
        density : float
            Density of the surrounding air; kg/m³
        thrust : float
            Thrust; N
        
        Returns
        -------
        float
            Blade loading.

        """
        # Blade loading [-]
        return (
            thrust / 
            (density * self.tip_velocity ** 2 * self.number_of_blades 
            * self.chord * self.radius))


    def get_induced_velocity(self, density, v_inf, alpha, thrust):
        """
        Calculate the induced velocity iteratively in forward flight (not valid 
        for low sink rate in axial flight). [p.253]

        Parameters
        ----------
        density : float
            Density of the surrounding air; kg/m³
        v_inf : float
            Inflow velocity; m/s
        alpha : float
            Angle of attack; rad
        thrust : float
            Thrust; N
        
        Returns
        -------
        float
            Induced velocity; m/s

        """
        area = np.pi * self.radius ** 2
        
        def induced_velocity(v_i):
            U = np.sqrt(
                v_inf ** 2 - 2 * v_inf * v_i * np.sin(alpha) + v_i ** 2)
            return thrust / (2 * density * area * U)

        # Induced velocity [m/s]
        return optimize.fixed_point(induced_velocity, 1)


    def get_profile_power(self, density, advance_ratio):
        """
        Calculate the profile power. [p.258]
        
        Parameters
        ----------
        density : float
            Density of the surrounding air; kg/m³
        advance_ratio : float
            Advance ratio.
        
        Returns
        -------
        float
            Profile power; W

        """
        # Profile power [W]
        return (
            1 / 8 * density * self.tip_velocity ** 3 * self.solidity 
            * self.zero_lift_drag_coeff * np.pi * self.radius ** 2
            * (1 + 4.65 * advance_ratio ** 2))


    def get_climb_power(self, flight_speed, gamma, weight):
        """
        Calculate the climb power due to the change in potential energy.
        
        Parameters
        ----------
        flight_speed : float
            Flight speed; m/s
        gamma : float
            Climb angle; rad
        weight : float
            Aircraft weight; N
        
        Returns
        -------
        float
            Climb power; W

        """
        # Climb power [W]
        return weight * flight_speed * np.sin(gamma)


    def get_min_power_radius(self, density, thrust):
        """
        Calculate the optimal rotor radius with respect to induced and profile 
        power in hover. [p.103]
        
        Parameters
        ----------
        density : float
            Density of the surrounding air; kg/m³
        thrust : float
            Thrust; N
        
        Returns
        -------
        float
            Optimal radius; m
        
        """
        # Optimal radius [m]
        return (
            1 / self.tip_velocity 
            * np.sqrt(2 * thrust / (density * np.pi)) 
            * (self.kappa / (self.solidity * self.zero_lift_drag_coeff)) 
            ** (1 / 3))


    def get_figure_of_merit(self, ideal_induced_power, profile_power):
        """
        Calculate the figure of merit.
        
        Parameters
        ----------
        ideal_induced_power : float
            Ideal induced power; W
        profile_power : float
            Profile power; W

        Returns
        -------
        float
            Figure of merit

        """
        # Figure of merit [-]
        return (
            ideal_induced_power / 
            (self.kappa * ideal_induced_power + profile_power))


    def in_ground_effect(self, induced_power, height):
        """ 
        Calculate the induced power in ground effect (IGE) according to Hayden. 
        [p.288]

        Parameters
        ----------
        induced_power : float
            Induced power, W
        height : float
            Height of the landing gear above ground; m
        
        Returns
        -------
        float
            Induced power in ground effect; W
        
        """
        z_R = self.installation_height + height
        k_G = 1 / (0.9926 + 0.0379 * (2 * self.radius / z_R) ** 2)

        # Induced power in ground effect [W]
        return k_G * induced_power
    