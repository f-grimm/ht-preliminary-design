# -*- coding: utf-8 -*-
"""
Created on 2022-02-12

@author: Fabian Grimm (f.grimm@tum.de)
"""

from aircraft import Aircraft
from mission import Mission
from rotor import Rotor
import yaml
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

class Helicopter(Aircraft):
    """ Aircraft sub-class for conventional helicopter configurations.

    Requires:
        Main rotor
        Tail rotor
    """
    def __init__(self, filename: str):      
        # Initialize aircraft level attributes
        Aircraft.__init__(self, filename)

        # Load data from YAML file
        with open('data/configurations/' + filename + '.yaml') as file:
            data = yaml.safe_load(file)

        # Add main and tail rotor
        self.main_rotor = Rotor(data['Main rotor'])
        self.tail_rotor = Rotor(data['Tail rotor'])


    def preliminary_design(self, mission: Mission, segment: int, logs=True):
        """
        """
        # Set the mission segment as defined in _main
        mission.set_mission_segment(segment)
        
        # Estimate MTOW
        self.initial_mtow_estimation(mission)

        # Initialize design loop
        mtow_list = [self.mtow]
        counter = 0
        error = 1

        while error > 1e-4:

            # Sizing
            self.main_rotor_sizing(mission)
            self.tail_rotor_sizing()
            powers, flight_state = self.refined_performance(mission)
            masses, empty_weight_comp = self.mass_estimation(powers, mission)
            
            # Check for convergence
            mtow_list.append(self.mtow)
            error = abs((mtow_list[-1] - mtow_list[-2]) / mtow_list[-1])
            counter += 1
            if counter > 1e3: raise ValueError('No convergence.')

        if logs:
            print(
                '\nPreliminary Design\n' + '-' * 18 + '\n\n'
                + f'Aicraft: {self.name}\n'
                + f'Mission: {mission.name}\n'
                + f'Segment: {mission.segment} ' 
                    + f'({mission.height:.0f}m, {mission.duration:.2f}h, '
                    + f'{mission.flight_speed:.0f}m/s)\n\n'
                + f'Iterations: {len(mtow_list) - 1}\n\n'
                + f' - Main rotor radius: {self.main_rotor.radius:14.2f} m\n'
                + f' - Tail rotor radius: {self.tail_rotor.radius:14.2f} m\n'
                + f' - Drag: {flight_state["drag"]:27.2f} N\n'
                + f' - Thrust: {flight_state["thrust"]:25.2f} N\n'
                + f' - Angle of attack: '
                    + f'{flight_state["alpha"] * 180 / np.pi:16.2f} deg\n'
                + f' - Advance ratio: '
                    + f'{flight_state["advance ratio"]:18.2f} -\n'
                + f' - Induced velocity: '
                    + f'{flight_state["induced velocity"]:15.2f} m/s\n'
                + f' - Total power: {(powers["total"] * 1e-3):20.2f} kW\n'
                + f' - Empty weight: {masses["empty weight"]:19.2f} kg\n'
                + f' - MTOW: {self.mtow:27.2f} kg\n')

            # Plots
            mission.plot_mission()
            self.plot_powers(mission.density, 80)
            self.plot_mtow_convergence(mtow_list)
            self.plot_masses(masses)
            self.plot_empty_weight_comp(empty_weight_comp)
            plt.show()


    def initial_mtow_estimation(self, mission: Mission):
        """ Estimate the initial maximum take-off weight based on the mission 
        profile; (SFC and empty weight ratio are constant). [pp.90-91]
        """
        fuel_mass = (self.sfc * self.power_available * mission.duration)
        useful_load = (fuel_mass + mission.payload + mission.crew_mass)

        # MTOW [kg]
        self.mtow = useful_load / (1 - self.empty_weight_ratio)


    def main_rotor_sizing(self, mission: Mission):
        """ Determine the main rotor radius; optimized for min. power in hover 
        [pp.98-172]
        """
        weight = self.mtow * self.gravity
        
        # Main rotor radius [m]
        self.main_rotor.radius = self.main_rotor.get_min_power_radius(
            mission.density, thrust=weight)


    def tail_rotor_sizing(self):
        """ Determine the tail rotor radius according to Layton [pp.204-213]
        """
        # Tail rotor radius [m]
        self.tail_rotor.radius = 0.4 * np.sqrt(2.2 * self.mtow * 1e-3)


    def refined_performance(self, mission: Mission):
        """ Determine the power comsumption based on the current mission 
        segment [pp.251-321]
        """
        # Flight state
        drag = self.get_fuselage_drag(mission.density, mission.flight_speed)
        alpha = self.get_angle_of_attack(drag, mission.climb_angle)
        advance_ratio = mission.flight_speed / self.main_rotor.tip_velocity
        thrust = self.get_thrust(
            drag, alpha, mission.climb_angle, advance_ratio)
        induced_velocity = self.main_rotor.get_induced_velocity(
                mission.density, mission.flight_speed, alpha, thrust)

        # Power calculation
        induced_power = self.main_rotor.kappa * thrust * induced_velocity
        profile_power = self.main_rotor.get_profile_power(
            mission.density, advance_ratio)
        parasite_power = self.get_parasite_power(
            mission.density, mission.flight_speed)
        climb_power = self.get_climb_power(
            mission.flight_speed, mission.climb_angle)
        main_rotor_power = (induced_power + profile_power + parasite_power 
                            + climb_power)
        tail_rotor_power = self.tail_rotor.power_fraction * main_rotor_power
        transmission_losses = (((1 / self.eta_transmission) - 1) 
                               * main_rotor_power)
        total_power = (main_rotor_power + tail_rotor_power 
                       + transmission_losses + self.accessory_power)

        # Engine requirement [p.298]
        temperature_ratio = mission.temperature / 288.15
        pressure_ratio = mission.pressure / 101325
        power_msl = (total_power * np.sqrt(temperature_ratio) 
                     / pressure_ratio)

        # Powers [W] and flight state [-, N, rad, N, m/s]
        return(
            {'total': total_power, 'msl': power_msl, 'induced': induced_power,
             'profile': profile_power, 'tail rotor': tail_rotor_power,  
             'climb': climb_power, 'transmission': transmission_losses, 
             'parasite': parasite_power, 'accessory': self.accessory_power},
            {'advance ratio': advance_ratio, 'drag': drag, 'alpha': alpha, 
             'thrust': thrust, 'induced velocity': induced_velocity})


    def mass_estimation(self, powers: dict, mission: Mission):
        """ Determine the fuel mass, empty weight, and new maximum take-off 
        weight based on the required power [pp.324-325]
        """
        fuel_mass = self.sfc * 1.1 * powers['total'] * mission.duration
        empty_weight_comp = self.empty_weight_estimation(
            fuel_mass, power=powers['msl'])
        empty_weight = sum([m for m in empty_weight_comp.values() if m > 0])
        
        # MTOW [kg]
        self.mtow = (empty_weight + fuel_mass + mission.payload 
                     + mission.crew_mass)

        # Masses and empty weight components [kg]
        return (
            {'empty weight': empty_weight, 'fuel': fuel_mass, 
             'payload': mission.payload, 'crew': mission.crew_mass},
            empty_weight_comp)


    def empty_weight_estimation(self, fuel_mass, power):
        """ Estimate empty weight components of medium helicopters. 
        [pp.408-416]
        """
        # Wetted fuselage surface [m^2] and chord [m]
        wetted_surface = 59.09386 * np.exp(0.0000194463 * self.mtow)
        chord = self.main_rotor.get_chord()
        
        # Mass estimation [kg]
        m_mr = (33 * self.main_rotor.radius * chord 
                * self.main_rotor.number_of_blades + 16)
        m_tr = 0.003942 * self.mtow + 5.66
        m_fus_t = 0.11907 * self.mtow - 66.666
        m_prop = 1.83 * (133.8 + 0.1156 * power * 1e-3)
        m_drive = (0.00000166 * (0.9 * self.mtow) ** 2 
                   + 0.087780096 * self.mtow - 113.81241656)
        m_tanks = 164.751 * np.log(fuel_mass / 2.948) - 751.33
        m_fcs = 95.6368 * np.exp(0.000111114 * self.mtow)
        m_i = 25.444 * np.log(power / 0.7457e3) - 141.62
        m_hyd = 0.003258 * self.mtow + 5.24
        m_el = 218.496 * np.log(wetted_surface / 0.092903) - 1267.49
        m_av = 113.4 + self.special_equipment
        m_furn = 0.854 * wetted_surface + 9.98 * self.number_of_seats - 4.54
        m_ac_ai = 55.542 * np.log(10.7369 * wetted_surface) - 331.21
        m_l_h = 38
        
        # Landing gear
        if self.landing_gear_type == 'Skids':
            m_lg = (0.011113004 * (0.9 * self.mtow / 0.453592) ** 0.8606 
                    * self.main_rotor.number_of_blades ** 0.8046)
        elif self.landing_gear_type == 'Rigid wheels':
            m_lg = (0.187333496 * (0.9 * self.mtow / 0.453592) ** 0.6662 
                    * self.number_of_legs ** 0.536)
        elif self.landing_gear_type == 'Retractable wheels':
            m_lg = (0.187333496 * (0.9 * self.mtow / 0.453592) ** 0.6662 
                    * 2 ** 0.1198 * self.number_of_legs ** 0.536)
        else: raise ValueError('Invalid landing gear type.')

        # Empty weight components [kg]
        return {
            'main rotor': m_mr, 'tail rotor': m_tr, 'fuselage, tail': m_fus_t, 
            'landing gear': m_lg, 'engines': m_prop, 'transmission': m_drive, 
            'fuel tanks': m_tanks, 'flight control systems': m_fcs, 
            'instruments': m_i, 'hydraulic systems': m_hyd, 
            'electrical systems': m_el, 'avionics': m_av, 'furnishing': m_furn, 
            'air conditioning, anti-ice': m_ac_ai, 'loading, handling': m_l_h}


    def plot_powers(self, density, max_velocity):
        """ Plot powers over a range of horizontal flight speeds for the 
        current MTOW. Level flight approximation is applied.
        """
        # List of speeds (x-axis)
        speeds = np.linspace(0, max_velocity, 100)

        # Initialize parameters
        P_i, P_0, P_p, P_tr, P_tl, P_a, P = [], [], [], [], [], [], []
        weight = self.mtow * self.gravity

        for V in speeds:

            # Flight state
            advance_ratio = V / self.main_rotor.tip_velocity
            induced_velocity = self.main_rotor.get_induced_velocity_level(
                density, V, thrust=weight)

            # Power calculation [kW]
            P_i.append(self.main_rotor.kappa * weight * induced_velocity 
                       * 1e-3)
            P_0.append(self.main_rotor.get_profile_power(
                density, advance_ratio) * 1e-3)
            P_p.append(self.get_parasite_power(density, V) * 1e-3)
            main_rotor_power = P_i[-1] + P_0[-1] + P_p[-1]
            P_tr.append(self.tail_rotor.power_fraction * main_rotor_power)
            P_tl.append(((1 / self.eta_transmission) - 1) * main_rotor_power)
            P_a.append(self.accessory_power * 1e-3)
            P.append(main_rotor_power + P_tr[-1] + P_tl[-1] + P_a[-1])

        # Plot
        fig, ax = plt.subplots(figsize=(8,5))
        ax.set_title('Power in level flight', fontweight='bold')
        ax.plot(speeds, P, label=r'$P$ Total power')
        ax.plot(speeds, P_i, label=r'$P_i$ Induced power')
        ax.plot(speeds, P_0, label=r'$P_0$ Profile power')
        ax.plot(speeds, P_p, label=r'$P_p$ Parasite power')
        ax.plot(speeds, P_tr, label=r'$P_{TR}$ Tail rotor power')
        ax.plot(speeds, P_tl, label=r'$P_{tl}$ Transmission losses')
        ax.plot(speeds, P_a, label=r'$P_a$ Accessory power')
        ax.set_xlabel('Velocity [m/s]')
        ax.set_ylabel('Power [kW]')
        fig.tight_layout()
        ax.legend()
        ax.grid()


    def plot_mtow_convergence(self, mtow_list: list):
        """ Plot MTOW over the iterations within the design loop.
        """
        # Figure, plot
        fig, ax = plt.subplots(figsize=(8,5))
        ax.set_title('MTOW convergence', fontweight='bold')
        ax.plot(range(len(mtow_list)), mtow_list, 'k')
        ax.set_xlabel('Iterations')
        ax.set_ylabel('MTOW [kg]')
        fig.tight_layout()
        ax.grid()


    def plot_masses(self, masses: dict):
        """ Plot mass composition (empty weight, fuel, crew, payload) in a pie 
        chart.
        """
        # Data
        labels = [name.capitalize() for name in masses.keys()]
        sizes = masses.values()
        explode = (0, 0, 0.2, 0)

        # Figure, plot
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.set_title('Maximum take-off weight', fontweight='bold', pad=15)
        ax.set_prop_cycle('color', plt.cm.get_cmap('Set2').colors)
        ax.pie(
            sizes, explode=explode, labels=labels, autopct='%1.1f%%',
            shadow=False, startangle=90, wedgeprops=dict(width=0.4),
            pctdistance=0.80)
        ax.margins(0.2, 0.2)
        ax.axis('equal')


    def plot_empty_weight_comp(self, empty_weight_comp: dict):
        """ Plot empty weight components in a pie chart.
        """
        # Data
        labels = list(empty_weight_comp.keys())
        sizes = [m if m > 0 else 0 for m in empty_weight_comp.values()]
        total = sum(sizes)
        label_first_n = 6
        for i in range(len(labels)):
            labels[i] = (labels[i].capitalize() 
                        + f' ({sizes[i] / total * 100:.1f}%)')

        # Figure, plot
        fig = plt.figure(figsize=(9.5, 4))
        ax = fig.add_subplot(121)
        ax.set_title('Empty weight', fontweight='bold', pad=15)
        ax.set_prop_cycle(
            'color', plt.cm.get_cmap('Set3').colors 
            + plt.cm.get_cmap('Pastel1').colors)
        pie = ax.pie(
            sizes, labels=labels[:label_first_n] + [''] 
            * (len(labels) - label_first_n), autopct=None, shadow=False, 
            startangle=90, wedgeprops=dict(width=0.4), pctdistance=0.80)
        ax2 = fig.add_subplot(122)
        ax2.axis('off')
        ax2.legend(
            pie[0][:(label_first_n - 1):-1], labels[:(label_first_n - 1):-1], 
            loc='center right')
        plt.subplots_adjust(wspace=-0.2)
        ax.margins(0.2, 0.2)
        ax.axis('equal')


