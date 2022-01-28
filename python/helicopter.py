
from aircraft import Aircraft
from rotor import Rotor
from engines import Engines
from fuselage import Fuselage
from landing_gear import LandingGear
from mission import Mission, atmosphere
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm

class Helicopter(Aircraft):
    """
    Aircraft sub-class for conventional helicopter configurations containing 
    methods for preliminary design and analysis. 
    
    Note
    ----
    Requires *Main rotor*, *Tail rotor*, *Engines*, *Fuselage*, and 
    *Landing Gear* section in the configuration file.
    
    Attributes
    ----------
    main_rotor : Rotor
        Main rotor.
    tail_rotor : Rotor
        Tail rotor.
    engines : Engines
        Turboshaft engines.
    fuselage : Fuselage
        Fuselage.
    landing_gear : LandingGear
        Landing Gear.
    
    """

    def __init__(self, filename: str):      
        """
        Load the configuration. 
        
        Parameters
        ----------
        filename : str
            Path and name of the YAML file containing the configuration.
        
        """
        # Initialize aircraft level attributes
        Aircraft.__init__(self, filename)

        # Load data from .yaml file
        with open(filename) as file:
            data = yaml.safe_load(file)

        # Add components
        self.main_rotor   = Rotor(data['Main rotor'])
        self.tail_rotor   = Rotor(data['Tail rotor'])
        self.engines      = Engines(data['Engines'])
        self.fuselage     = Fuselage(data['Fuselage'])
        self.landing_gear = LandingGear(data['Landing gear'])


    def preliminary_design(
            self, mission: Mission, radius_range: tuple, chord_range: tuple, 
            conv_tol=1e-4, use_ew_models=False, status=True):
        """
        Generate a dataframe for given ranges of the rotor radius and chord 
        length, containing maximum take-off weight, fuel consumption, and 
        other performance parameters.

        Parameters
        ----------
        mission : Mission
            Mission profile.
        radius_range : tuple(min, max, N)
            Range of the rotor radius variation; m
        chord_range : tuple(min, max, M)
            Range of the chord length variation; m
        conv_tol : float, optional
            Convergence tolerance of the sizing loop.
        use_ew_models : bool, optional
            Use empirical empty weight models, otherwise the ratio to the 
            maximum take-off weight is assumed constant.
        status : bool, optional
            Print status updates, default True.
        
        Returns
        -------
        pd.DataFrame
            Dataframe containing the results of :meth:`sizing_loop` as rows.

        """
        # Initialize empty dataframe
        df = pd.DataFrame()

        # Design space
        for radius in np.linspace(*radius_range):
            for chord in np.linspace(*chord_range):

                # Set rotor parameters
                self.main_rotor.radius = radius
                self.main_rotor.chord = chord

                # Intitial MTOW
                self.mtow = self.initial_mtow_estimation(
                    sum(mission.duration), min(mission.payload), 
                    min(mission.crew_mass))

                # Converged MTOW and fuel consumption
                result = self.sizing_loop(mission, conv_tol, use_ew_models)
                
                # Append to dataframe as new row
                new_row = pd.DataFrame([result])
                df = pd.concat([df, new_row], ignore_index=True)

                # Print status
                if status: print(
                    'Calculating... '
                    f'{len(df)}/{radius_range[2] * chord_range[2]}')

        if status: print('\nDone!\n')

        # Return all results
        return df
        

    def sizing_loop(self, mission: Mission, conv_tol, use_ew_models):
        """
        Determine the maximum take-off weight and performance iteratively for a 
        given mission.

        Parameters
        ----------
        mission : Mission
            Mission profile.
        conv_tol : float
            Convergence tolerance.
        use_ew_models : bool
            Use empirical empty weight models, otherwise the ratio to the 
            maximum take-off weight is assumed constant.
        
        Returns
        -------
        dict
            Dictionary containing the radius, chord length, maximum take-off
            weight, fuel mass, empty weight ratio, solidity, disc loading, and
            blade loading.

        """        
        # Initialize sizing loop
        converged = False
        
        while not converged:
            
            # Perform mission with the current MTOW
            values_per_segment = self.evaluate(mission)

            # New mass estimation
            fuel_total = (
                sum(values_per_segment['fuel']) * self.engines.fuel_safety)
            max_power_msl = max([
                performance['power msl'] for performance in 
                values_per_segment['performance']])
            masses = self.mass_estimation(
                max_power_msl, fuel_total, mission.crew_mass[0], 
                mission.payload[0], use_ew_models)
            new_mtow = sum(masses.values())
            
            # Check for convergence
            converged = abs((new_mtow - self.mtow) / new_mtow) < conv_tol
            
            # Update MTOW
            self.mtow = new_mtow

        # Final variables
        empty_weight_ratio = masses['empty weight'] / new_mtow
        max_thrust = max([
            fs['thrust'] for fs in values_per_segment['flight state']])
        max_disc_loading = self.main_rotor.get_disc_loading(max_thrust)
        max_blade_loading = max([
            fs['blade loading'] for fs in values_per_segment['flight state']])

        # Result
        return {
            'Radius [m]': self.main_rotor.radius, 
            'Chord [m]': self.main_rotor.chord, 'MTOW [kg]': new_mtow, 
            'Fuel [kg]': fuel_total, 'EW ratio': empty_weight_ratio, 
            'Solidity': self.main_rotor.solidity, 
            'DL [N/m²]': max_disc_loading, 'Blade ldg.': max_blade_loading}


    def evaluate(self, mission: Mission):
        """
        Evaluate the flight state, powers, and performance throughout the 
        mission for a fixed configuration. For each segment, the Gross Weight 
        at the beginning of it is assumed.
        
        Parameters
        ----------
        mission: Mission
            Mission profile.

        Returns
        -------
        dict
            Dictionary containing lists of the flight state, power composition, 
            performance, fuel mass, and gross weight throughout the mission.

        """
        # Initialize
        gross_weight_per_segment = []
        flight_state_per_segment = []
        powers_per_segment = []
        performance_per_segment = []
        fuel_mass_per_segment = []

        # Go through the mission
        for i in range(len(mission.duration)):

            # Determine the current aircraft mass
            gross_weight_per_segment.append(
                self.mtow - sum(fuel_mass_per_segment[:i]) 
                + mission.payload[i] - mission.payload[0] 
                + mission.crew_mass[i] - mission.crew_mass[0])

            # Flight state
            flight_state = self.flight_state(
                gross_weight_per_segment[i], mission.density[i], 
                mission.flight_speed[i], mission.climb_angle[i], 
                mission.gravity[i])
            flight_state_per_segment.append(flight_state)

            # Power composition
            powers = self.powers(flight_state)
            powers_per_segment.append(powers)

            # Performance and fuel calculation
            performance = self.performance(
                powers, mission.temperature[i], mission.pressure[i])
            performance_per_segment.append(performance)
            fuel_mass_per_segment.append(
                performance['fuel mass flow'] * mission.duration[i])

        # Values per segment
        return {
            'flight state': flight_state_per_segment, #: list[dict]
            'powers': powers_per_segment, #: list[dict]
            'performance': performance_per_segment, #: list[dict]
            'fuel': fuel_mass_per_segment, #: list[float]
            'GW': gross_weight_per_segment} #: list[float]


    def initial_mtow_estimation(self, duration, payload, crew_mass):
        """
        Estimate the initial maximum take-off weight based on the installed 
        power and total duration. Specific fuel consumption and empty weight 
        ratio are assumed constant. [pp.90-91]
        
        Parameters
        ----------
        duration : float
            Mission duration; h
        payload : float
            Payload; kg
        crew_mass : float
            Crew mass; kg
        
        """
        fuel_mass = (
            self.engines.sfc * self.engines.power_available 
            * self.engines.number_of_engines * duration)
        useful_load = fuel_mass + payload + crew_mass

        # MTOW [kg]
        return useful_load / (1 - self.empty_weight_ratio)


    def flight_state(self, gross_weight, density, velocity, gamma, gravity):
        """
        Determine flight state variables.
        
        Parameters
        ----------
        gross_weight : float
            Aircraft mass; kg
        density : float
            Density of the surrounding air; kg/m³
        velocity : float
            Flight speed; m/s
        gamma : float
            Climb angle; rad
        gravity : float
            Gravitational acceleration; m/s²
        
        Returns
        -------
        dict
            Dictionary containing the density, velocity, climb angle, advance 
            ratio, drag, angle of attack, thrust, induced velocity, download 
            factor, blade loading, gravity, and gross weight.
        
        """
        drag = self.fuselage.get_fuselage_drag(density, velocity)
        thrust_isolated, alpha = self.get_thrust_and_alpha(
            gross_weight, drag, gamma, gravity)
        advance_ratio = velocity / self.main_rotor.tip_velocity
        k_DL = self.fuselage.get_download_factor_in_flight(advance_ratio)
        thrust = thrust_isolated / (1 - k_DL)
        blade_loading = self.main_rotor.get_blade_loading(density, thrust)
        induced_velocity = self.main_rotor.get_induced_velocity(
            density, velocity, alpha, thrust)

        # Flight state
        return {
            'density': density, 'velocity': velocity, 'gamma': gamma, 
            'advance ratio': advance_ratio, 'drag': drag, 'alpha': alpha, 
            'thrust': thrust, 'induced velocity': induced_velocity,
            'download factor': k_DL, 'blade loading': blade_loading, 
            'gravity': gravity, 'GW': gross_weight}
    

    def powers(self, flight_state: dict):
        """
        Determine the power composition in the current flight state.
        
        Parameters
        ----------
        flight_state : dict
            Output of :meth:`flight_state`.
        
        Returns
        -------
        dict
            Dictionary containing the induced, profile, climb, parasite, tail 
            rotor, and accessory power as well as transmission losses. The sum 
            represents the total power.

        """
        density = flight_state['density']
        velocity = flight_state['velocity']
        gamma = flight_state['gamma']
        thrust = flight_state['thrust']
        induced_velocity = flight_state['induced velocity']
        advance_ratio = flight_state['advance ratio']
        weight = flight_state['GW'] * flight_state['gravity']

        # Power calculation        
        induced_power = self.main_rotor.kappa * thrust * induced_velocity
        profile_power = self.main_rotor.get_profile_power(
            density, advance_ratio)
        parasite_power = self.fuselage.get_parasite_power(density, velocity)
        climb_power = self.main_rotor.get_climb_power(velocity, gamma, weight)
        main_rotor_power = (
            induced_power + profile_power + parasite_power + climb_power)
        tail_rotor_power = self.tail_rotor.power_fraction * main_rotor_power
        transmission_losses = (
            ((1 / self.eta_transmission) - 1) * main_rotor_power)

        # Powers [W]
        return{
            'induced': induced_power, 'profile': profile_power, 
            'climb': climb_power, 'parasite': parasite_power, 
            'tail rotor': tail_rotor_power, 
            'transmission': transmission_losses, 
            'accessory': self.accessory_power}


    def performance(self, powers: dict, temperature, pressure):
        """
        Determine performance parameters aside from the individual power
        components.

        Parameters
        ----------
        powers : dict
            Output of :meth:`powers`.
        temperature : float
            Temperature of the surrounding air; K
        pressure : float
            Pressure of the surrounding air; Pa
        
        Returns
        -------
        dict
            Dictionary containing the total power, the power requirement at 
            mean sea level, specific fuel consumption, and the fuel mass flow.
        
        """
        # Total power
        total_power = sum(powers.values())

        # Engine requirement at mean sea level [p.298]
        temperature_ratio = temperature / 288.15
        pressure_ratio = pressure / 101325
        power_msl = total_power * np.sqrt(temperature_ratio) / pressure_ratio

        # Use power-dependent SFC if engine parameters A, B are provided
        if self.engines.a is not None and self.engines.b is not None:
            sfc = self.engines.get_sfc(
                temperature_ratio, pressure_ratio, total_power)
        else:
            sfc = self.engines.sfc

        # Fuel mass flow [kg/h]
        fuel_mass_flow = sfc * total_power

        # Performance
        return {
            'total power': total_power, 'power msl': power_msl, 'SFC': sfc,
            'fuel mass flow': fuel_mass_flow}


    def mass_estimation(
            self, power, fuel_mass, crew_mass, payload, use_ew_models):
        """
        Determine the maximum take-off weight composition.

        Parameters
        ----------
        power : float
            Power; W
        fuel_mass : float
            Fuel mass; kg
        crew_mass : float
            Crew mass; kg
        payload : float
            Payload; kg
        use_ew_models : bool
            Use empirical empty weight models, otherwise the ratio to the 
            maximum take-off weight is assumed constant.

        Returns
        -------
        dict
            Dictionary containing the empty weight, fuel, payload, and crew 
            mass. The sum represents the maximum take-off weight.

        """
        if use_ew_models:
            empty_weight_parts = self.empty_weight_estimation(fuel_mass, power)
            empty_weight = sum(
                [m for m in empty_weight_parts.values() if m > 0])
        else:
            empty_weight = (
                self.empty_weight_ratio * (fuel_mass + crew_mass + payload) / 
                (1 - self.empty_weight_ratio))
        
        # Masses [kg]
        return {
            'empty weight': empty_weight, 'fuel': fuel_mass, 
            'payload': payload, 'crew': crew_mass}


    def empty_weight_estimation(self, fuel_mass, power):
        """
        Estimate the empty weight components of medium helicopters. 
        [pp.408-416]
        
        Parameters
        ----------
        fuel_mass : float
            Fuel mass; kg, used for fuel tank sizing.
        power : float
            Total power; W, used for propulsion system and instrument sizing.
        
        Returns
        -------
        dict
            Dictionary containing the empty weight components main rotor, tail 
            rotor, fuselage and tail, landing gear, transmission, engines, fuel 
            tanks, flight control systems, hydraulic systems, avionics, 
            instruments, furnishing, air conditioning and anti-ice, loading and 
            handling, and electrical systems. The sum represents the total 
            empty weight.

        """
        # Wetted fuselage surface [m²]
        wetted_surface = 59.09386 * np.exp(0.0000194463 * self.mtow)
        
        # Mass estimation [kg]
        m_main_rotor = (
            33 * self.main_rotor.radius * self.main_rotor.chord 
            * self.main_rotor.number_of_blades + 16)
        m_tail_rotor = 0.003942 * self.mtow + 5.66
        m_fuselage_tail = 0.11907 * self.mtow - 66.666
        m_propulsion = 1.83 * (133.8 + 0.1156 * power * 1e-3)
        m_drive_sys = (
            0.00000166 * (0.9 * self.mtow) ** 2 
            + 0.087780096 * self.mtow - 113.81241656)
        m_fuel_tanks = 164.751 * np.log(fuel_mass / 2.948) - 751.33
        m_flight_control_sys = 95.6368 * np.exp(0.000111114 * self.mtow)
        m_instruments = 25.444 * np.log(power / 0.7457e3) - 141.62
        m_hydraulics = 0.003258 * self.mtow + 5.24
        m_electronics = 218.496 * np.log(wetted_surface / 0.092903) - 1267.49
        m_avionics = 113.4 + self.special_equipment
        m_furnishing = (
            0.854 * wetted_surface + 9.98 * self.fuselage.number_of_seats 
            - 4.54)
        m_ac_anti_ice = 55.542 * np.log(10.7369 * wetted_surface) - 331.21
        m_loading_handling = 38
        
        # Landing gear
        if self.landing_gear.type_ == 'Skids':
            m_landing_gear = (
                0.011113004 * (0.9 * self.mtow / 0.453592) ** 0.8606 
                * self.main_rotor.number_of_blades ** 0.8046)
        elif self.landing_gear.type_ == 'Rigid wheels':
            m_landing_gear = (
                0.187333496 * (0.9 * self.mtow / 0.453592) ** 0.6662 
                * self.landing_gear.number_of_legs ** 0.536)
        elif self.landing_gear.type_ == 'Retractable wheels':
            m_landing_gear = (
                0.187333496 * (0.9 * self.mtow / 0.453592) ** 0.6662 
                * 2 ** 0.1198 * self.number_of_legs ** 0.536)
        else: raise ValueError('Invalid landing gear type.')

        # Empty weight components [kg]
        return {
            'main rotor': m_main_rotor, 'tail rotor': m_tail_rotor, 
            'fuselage, tail': m_fuselage_tail, 'landing gear': m_landing_gear, 
            'transmission': m_drive_sys, 'engines': m_propulsion, 
            'fuel tanks': m_fuel_tanks, 
            'flight control systems': m_flight_control_sys, 
            'hydraulic systems': m_hydraulics, 'avionics': m_avionics, 
            'instruments': m_instruments, 'furnishing': m_furnishing, 
            'air conditioning, anti-ice': m_ac_anti_ice, 
            'loading, handling': m_loading_handling, 
            'electrical systems': m_electronics}


    def plot_mtow_convergence(self, mtow_list: list):
        """
        Plot the maximum take-off weight over the iterations.
        
        Parameters
        ----------
        mtow_list : list[float]
            List of the maximum take-off weight; kg
        
        """
        # Plot
        fig, ax = plt.subplots(figsize=(8,5))
        ax.set_title('MTOW convergence', fontweight='bold')
        ax.plot(range(len(mtow_list)), mtow_list, 'k')
        ax.set_xlabel('Iterations')
        ax.set_ylabel('MTOW [kg]')
        fig.tight_layout()
        ax.grid()


    def plot_power_curves(self, gross_weight, density, max_velocity=None):
        """
        Plot the individual power components over a range of horizontal flight 
        speeds.
        
        Parameters
        ----------
        gross_weight : float
            Aircraft mass; kg
        density : float
            Density of the surrounding air; kg/m³
        max_velocity : float, optional
            Upper limit for the forward flight speed to be considered in the 
            plot. If not provided, the limit is defined by the installed power.
        
        """
        # Initialize parameters
        P_i, P_0, P_p, P_tr, P_tl, P_a, P = [], [], [], [], [], [], []
        speeds = []
        max_reached = False
        V = 0

        while not max_reached:

            # Power calculation
            flight_state = self.flight_state(
                gross_weight, density, V, gamma=0, gravity=9.81)
            powers = self.powers(flight_state)
            P_i.append(powers['induced'] * 1e-3)
            P_0.append(powers['profile'] * 1e-3)
            P_p.append(powers['parasite'] * 1e-3)
            P_tr.append(powers['tail rotor'] * 1e-3)
            P_tl.append(powers['transmission'] * 1e-3)
            P_a.append(powers['accessory'] * 1e-3)
            P.append(sum(powers.values()) * 1e-3)

            # Increment velocity
            speeds.append(V)
            V += 1.0

            # Exit when the provided max. velocity is exceeded
            if max_velocity is not None:
                max_reached = V > max_velocity
            
            # Exit when the installed power is exceeded
            else: max_reached = (
                P[-1] > self.engines.power_available 
                        * self.engines.number_of_engines * 1e-3)

        if max_velocity is None and len(P) == 1:
            print('\nWarning: Installed power not sufficient for hover.')

        # Plot
        fig, ax = plt.subplots(figsize=(8,5))
        ax.set_title(f'Power in level flight', fontweight='bold')
        ax.plot(speeds, P, label=r'$P$, total power')
        ax.plot(speeds, P_i, label=r'$P_i$, induced power')
        ax.plot(speeds, P_0, label=r'$P_0$, profile power')
        ax.plot(speeds, P_p, label=r'$P_p$, parasite power')
        ax.plot(speeds, P_tr, label=r'$P_{TR}$, tail rotor power')
        ax.plot(speeds, P_tl, label=r'$P_{tl}$, transmission losses')
        ax.plot(speeds, P_a, label=r'$P_a$, accessory power')
        ax.set_xlabel('Velocity [m/s]')
        ax.set_ylabel('Power [kW]')
        fig.tight_layout()
        ax.legend()
        ax.grid()


    def plot_fuel_curve(
            self, gross_weight, density, temperature, pressure, 
            max_velocity=None):
        """
        Plot the fuel mass flow and total power over a range of horizontal 
        flight speeds. Annotations show the velocities for best range and 
        endurance.

        Parameters
        ----------
        gross_weight : float
            Aircraft mass; kg
        density : float
            Density of the surrounding air; kg/m³
        temperature : float
            Temperature of the surrounding air; K
        pressure : float
            Pressure of the surrounding air; Pa
        max_velocity : float, optional
            Upper limit for the forward flight speed to be considered in the 
            plot. If not provided, the limit is defined by the installed power.
        
        """
        # Initialize parameters
        P, m_dots, speeds = [], [], []
        max_reached = False
        V = 0

        while not max_reached:

            # Power and fuel calculation
            flight_state = self.flight_state(
                gross_weight, density, V, gamma=0, gravity=9.81)
            powers = self.powers(flight_state)
            performance = self.performance(powers, temperature, pressure)
            P.append(sum(powers.values()) * 1e-3)
            m_dots.append(performance['fuel mass flow'])

            # Increment velocity
            speeds.append(V)
            V += 1.0

            # Exit when the provided max. velocity is exceeded
            if max_velocity is not None:
                max_reached = V > max_velocity
            
            # Exit when the installed power is exceeded
            else: max_reached = (
                P[-1] > self.engines.power_available 
                        * self.engines.number_of_engines * 1e-3)

        # Plot
        fig, ax = plt.subplots(figsize=(8,5))
        ax2 = ax.twinx()
        ax.set_title('Fuel mass flow', fontweight='bold')
        ax.plot(speeds, P, color='tab:blue', label='Total power')
        ax.set_xlabel('Velocity [m/s]')
        ax.set_ylabel('Power [kW]')
        ax.legend(loc='lower left')
        ax.set_ylim(bottom=0)
        ax2.plot(speeds, m_dots, color='tab:brown', label='Fuel mass flow')
        ax2.set_ylabel('Mass flow [kg/h]')
        ax2.legend(loc='lower right')
        ax2.set_ylim(bottom=0, top=m_dots[-1] / 0.8)
        fig.tight_layout()
        ax.grid()

        # Best range and endurance
        m_dot_V = [m_dot / V for m_dot, V in zip(m_dots[1:], speeds[1:])]
        index_best_range = m_dot_V.index(min(m_dot_V)) + 1
        index_best_endurance = m_dots.index(min(m_dots))
        V_best_range = speeds[index_best_range]
        V_best_endurance = speeds[index_best_endurance]
        ax2.annotate(
            r'$V_{best\,range}$ = ' + f'{V_best_range:.0f} m/s', 
            xy=(V_best_range, m_dots[index_best_range]), 
            xytext=(0, -40), textcoords='offset points', ha='center',
            arrowprops=dict(arrowstyle='->', color='black'))
        ax2.annotate(
            r'$V_{best\,endurance}$ = ' + f'{V_best_endurance:.0f} m/s', 
            xy=(V_best_endurance, m_dots[index_best_endurance]), 
            xytext=(0, -40), textcoords='offset points', ha='center',
            arrowprops=dict(arrowstyle='->', color='black'))
    

    def plot_power_sweep_gw(self, gw_range: tuple, density, max_velocity=75):
        """
        Plot power curves for a range of Gross Weights.

        Parameters
        ----------
        gw_range : tuple(min, max, N)
            Gross Weight range; kg
        density : float
            Density of the surrounding air; kg/m³
        max_velocity : float, optional
            Maximum velocity to be considered in the plot; m/s

        """
        # Initialize parameters and figure
        speeds = np.linspace(0, max_velocity, 100)
        gross_weights = np.linspace(*gw_range)
        fig, ax = plt.subplots(figsize=(8,5))
        cmap = cm.coolwarm(np.linspace(0, 1, gw_range[2]))

        for i, gross_weight in enumerate(gross_weights):

            # Reset power curve
            P = []
            
            for V in speeds:

                # Power calculation
                flight_state = self.flight_state(
                    gross_weight, density, V, gamma=0, gravity=9.81)
                powers = self.powers(flight_state)
                P.append(sum(powers.values()) * 1e-3)

            # Plot curve
            label = f'GW = {gross_weight:.0f} kg'
            ax.plot(speeds, P, color=cmap[i], label=label)

        # Plot settings
        ax.set_title(f'Influence of the Gross Weight', fontweight='bold')
        ax.set_xlabel('Velocity [m/s]')
        ax.set_ylabel('Power [kW]')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1])
        ax.set_ylim(bottom=0)
        fig.tight_layout()
        ax.grid()


    def plot_power_sweep_drag_area(
            self, drag_area_range: tuple, gross_weight, density, 
            max_velocity=75):
        """
        Plot power curves for a range of drag areas.

        Parameters
        ----------
        drag_area_range : tuple(min, max, N)
            Drag area range; kg
        gross_weight : float
            Aircraft mass; kg
        density : float
            Density of the surrounding air; kg/m³
        max_velocity : float, optional
            Maximum velocity to be considered in the plot; m/s

        """
        # Initialize parameters and figure
        restore = self.fuselage.drag_area
        speeds = np.linspace(0, max_velocity, 100)
        drag_areas = np.linspace(*drag_area_range)
        fig, ax = plt.subplots(figsize=(8,5))
        cmap = cm.coolwarm(np.linspace(0, 1, drag_area_range[2]))

        for i, drag_area in enumerate(drag_areas):

            # Reset power curve
            P = []

            # Parameter variation
            self.fuselage.drag_area = drag_area
            
            for V in speeds:

                # Power calculation
                flight_state = self.flight_state(
                    gross_weight, density, V, gamma=0, gravity=9.81)
                powers = self.powers(flight_state)
                P.append(sum(powers.values()) * 1e-3)

            # Plot curve
            label = f'Drag area = {drag_area:.2f} m²'
            ax.plot(speeds, P, color=cmap[i], label=label)

        # Restore
        self.fuselage.drag_area = restore

        # Plot settings
        ax.set_title(f'Influence of the drag area', fontweight='bold')
        ax.set_xlabel('Velocity [m/s]')
        ax.set_ylabel('Power [kW]')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1])
        ax.set_ylim(bottom=0)
        fig.tight_layout()
        ax.grid()


    def plot_power_sweep_height(
            self, height_range: tuple, gross_weight, max_velocity=75):
        """
        Plot power curves for a range of heights.

        Parameters
        ----------
        height_range : tuple(min, max, N)
            Height range; kg
        gross_weight : float
            Aircraft mass; kg
        max_velocity : float, optional
            Maximum velocity to be considered in the plot; m/s

        """
        # Initialize parameters and figure
        speeds = np.linspace(0, max_velocity, 100)
        heights = np.linspace(*height_range)
        fig, ax = plt.subplots(figsize=(8,5))
        cmap = cm.coolwarm(np.linspace(0, 1, height_range[2]))

        for i, height in enumerate(heights):

            # Reset power curve
            P = []

            # Parameter variation
            density, _, _ = atmosphere(height, temp_offset=0)
            
            for V in speeds:

                # Power calculation
                flight_state = self.flight_state(
                    gross_weight, density, V, gamma=0, gravity=9.81)
                powers = self.powers(flight_state)
                P.append(sum(powers.values()) * 1e-3)

            # Plot curve
            label = f'h = {height:.0f} m'
            ax.plot(speeds, P, color=cmap[i], label=label)

        # Plot settings
        ax.set_title(f'Influence of the height', fontweight='bold')
        ax.set_xlabel('Velocity [m/s]')
        ax.set_ylabel('Power [kW]')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1])
        ax.set_ylim(bottom=0)
        fig.tight_layout()
        ax.grid()


    def plot_powers_pie(self, powers: dict):
        """
        Plot the power composition as a pie chart.
        
        Parameters
        ----------
        powers : dict
            Output of :meth:`powers`.
        
        """
        # Data, consider only positive elements
        labels, sizes = [], []
        for name, P in powers.items():
            if P > 0: 
                labels.append(name.capitalize())
                sizes.append(P)
            else: 
                labels.append('')
                sizes.append(0)
                
        # Plot
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.set_title('Total power', fontweight='bold', pad=15)
        ax.set_prop_cycle('color', plt.cm.get_cmap('Set2').colors[1:])
        ax.pie(
            sizes, labels=labels, shadow=False, startangle=90, 
            wedgeprops=dict(width=0.4), pctdistance=0.80,
            autopct=lambda p:'{:.1f}%'.format(round(p)) if p > 0 else '')
        ax.margins(0.2, 0.2)
        ax.axis('equal')


    def plot_masses_pie(self, masses: dict):
        """
        Plot the maximum take-off weight composition as a pie chart.

        Parameters
        ----------
        masses : dict
            Output of :meth:`mass_estimation`.
        
        """
        # Data, consider only positive elements
        labels, sizes = [], []
        explode = (0, 0, 0.2, 0)
        for name, m in masses.items():
            if m > 0:
                labels.append(name.capitalize())
                sizes.append(m)
            else:
                labels.append('')
                sizes.append(0)

        # Plot
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.set_title('Maximum take-off weight', fontweight='bold', pad=15)
        ax.set_prop_cycle('color', plt.cm.get_cmap('Set2').colors)
        ax.pie(
            sizes, explode=explode, labels=labels, shadow=False, startangle=90, 
            wedgeprops=dict(width=0.4), pctdistance=0.80,
            autopct=lambda p:'{:.1f}%'.format(round(p)) if p > 0 else '')
        ax.margins(0.2, 0.2)
        ax.axis('equal')


    def plot_empty_weight_pie(self, empty_weight_parts: dict):
        """
        Plot the empty weight composition as a pie chart.

        Parameters
        ----------
        empty_weight_parts : dict
            Output of :meth:`empty_weight_estimation`.
        
        """
        # Data, consider only positive elements
        labels, sizes = [], []
        for name, m in empty_weight_parts.items():
            if m > 0:
                labels.append(name.capitalize())
                sizes.append(m)
            else:
                labels.append('')
                sizes.append(0)
        
        # Plot
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.set_title('Empty weight', fontweight='bold', pad=15)
        ax.set_prop_cycle(
            'color', plt.cm.get_cmap('Set3').colors 
            + plt.cm.get_cmap('Pastel1').colors)
        ax.pie(
            sizes, labels=labels, shadow=False, startangle=90, 
            wedgeprops=dict(width=0.4), pctdistance=0.80,
            autopct=lambda p:'{:.1f}%'.format(round(p)) if p > 0 else '')
        ax.margins(0.2, 0.2)
        ax.axis('equal')


    def plot_blade_loading(self, blade_loading, advance_ratio):
        """
        Plot the blade loading over the advance ratio, in order to evaluate 
        the stall margin. [p.142]

        Parameters
        ----------
        blade_loading : float
            Blade loading.
        advance_ratio : float
            Advance ratio.
        
        """
        # Empirical estimation for the rotor limit
        a = -0.16 / 0.6
        b = 0.16
        mu_limit = np.linspace(0.1, 0.5, 2)
        ct_sigma_limit = [a * mu + b for mu in mu_limit]

        if blade_loading > a * advance_ratio + b: print(
            '\nWarning: Blade loading exceeds the stall limit. '
            'Recommendation to increase the chord length.\n')

        # Plot
        fig, ax = plt.subplots(figsize=(6,4))
        ax.set_title('Blade loading', fontweight='bold')
        ax.plot(mu_limit, ct_sigma_limit, 'k-', label='Rotor limit')
        ax.plot(
            advance_ratio, blade_loading, 'o', color='tab:blue', 
            label=self.name, markersize=8)
        ax.set_xlabel(r'Advance ratio $\mu$ [-]')
        ax.set_ylabel(r'Blade loading  $C_T / \sigma$ [-]')
        fig.tight_layout()
        ax.legend()
        ax.grid()


    def plot_gross_weight_over_time(
            self, mission: Mission, fuel_mass_per_segment: list):
        """
        Plot the Gross Weight and fuel flow over the mission duration.

        Parameters
        ----------
        mission : Mission
            Mission profile.
        fuel_mass_per_segment : list[float]
            Fuel demand for each segment; kg
        
        """
        # Time, double elements for discontinuous lines
        time = (
            [0] + [t for t in np.cumsum(mission.duration) 
            for _ in (0, 1)][:-1])

        # Gross Weight
        gross_weight_list = [self.mtow]

        for i, fuel_mass in enumerate(fuel_mass_per_segment):
            gross_weight_list.append(gross_weight_list[-1] - fuel_mass)
            if i < len(mission.duration) - 1:
                gross_weight_list.append(
                    gross_weight_list[-1] + mission.payload[i + 1] 
                    - mission.payload[i] + mission.crew_mass[i + 1]
                    - mission.crew_mass[i]) 

        # Fuel flow
        fuel_flow = np.array(fuel_mass_per_segment) / mission.duration 
        fuel_flow_list = [m_dot for m_dot in fuel_flow for _ in (0, 1)]
        
        # Plot
        fig, ax = plt.subplots(figsize=(8,5))
        ax2 = ax.twinx()
        ax.set_title('Gross Weight', fontweight='bold')
        ax.plot(time, gross_weight_list, color='tab:green', label='GW')        
        ax.grid()
        ax.set_xlabel('Time [h]')
        ax.set_ylabel('Mass [kg]')
        ax.legend(loc='lower left')
        ax2.plot(time, fuel_flow_list, color='tab:brown', label='Fuel flow')
        ax2.set_ylabel('Mass flow [kg/h]')
        ax2.legend(loc='lower right')
        ax2.set_ylim(bottom=-10, top=max(fuel_flow) / 0.6)
        fig.tight_layout()


    def plot_results(self, df: pd.DataFrame, param: str):
        """
        Plot one of the preliminary design parameters as a 3D surface over
        rotor radius and chord length.

        Parameters
        ----------
        df : pd.DataFrame
            Output of :meth:`preliminary_design`.
        param : str
            z-axis parameter, key of :meth:`sizing_loop` output.

        """
        radius_column = df['Radius [m]'].to_numpy()
        chord_column = df['Chord [m]'].to_numpy()
        param_column = df[param].to_numpy()
        radius_list = np.unique(radius_column)
        chord_list = np.unique(chord_column)
        z = param_column.reshape((len(radius_list), len(chord_list)))
        x, y = np.meshgrid(chord_list, radius_list)
        
        # Plot 3D surface
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        surf = ax.plot_surface(
            x, y, z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        ax.set_title(param, fontweight='bold')
        ax.set_xlabel('Chord [m]')
        ax.set_ylabel('Radius [m]')
        fig.colorbar(surf, shrink=0.5, aspect=10, pad=0.1)
