# -*- coding: utf-8 -*-
"""
Created on 2022-01-30

@author: Fabian Grimm (f.grimm@tum.de)
"""

import numpy as np
import yaml
import rotor

class Aircraft:
    """
    """
    def __init__(self, filename: str):
        # Load data from YAML file
        with open('data/configurations/' + filename + '.yaml') as file:
            data = yaml.safe_load(file)

        self.name               = data['Name']
        self.power_available    = data['Engine'].get('Power available')
        self.accessory_power    = data['Engine'].get('Accessory power')
        self.sfc                = data['Engine'].get('SFC')
        self.download_factor    = data['Fuselage'].get('Download factor')
        self.drag_area          = data['Fuselage'].get('Drag area')
        self.landing_gear_type  = data['Landing gear'].get('Type')
        self.number_of_legs     = data['Landing gear'].get('Number of legs')
        self.empty_weight_ratio = data['Misc'].get('Empty weight ratio')
        self.eta_transmission   = data['Misc'].get('Transmission efficiency')
        self.special_equipment  = data['Misc'].get('Special equipment', 0)
        self.number_of_seats    = data['Misc'].get('Number of seats')
        self.gravity            = data['Misc'].get('Gravity', 9.81)
        

    def get_parasite_power(self, density, flight_speed):
        """ Calculate the parasite power created by the fuselage drag in 
        forward flight.
        """
        # Parasite power [W]
        return (0.5 * density * flight_speed ** 3 * self.drag_area)


    def get_climb_power(self, flight_speed, gamma):
        """ Claculate the climb power due to change in potential energy. Climb 
        angle gamma in rad.
        """
        weight = self.mtow * self.gravity

        # Climb power [W]
        return (weight * flight_speed * np.sin(gamma))


    def get_fuselage_drag(self, density, flight_speed):
        """ Calculate the fuselage drag.
        """
        # Fuselage drag [N]
        return (0.5 * density * flight_speed ** 2 * self.drag_area)


    def get_thrust(self, drag, alpha, gamma, advance_ratio):
        """ Calculate the required thrust in translatory motion. Angle of 
        attack alpha (relative to the flight path) and climb angle gamma in 
        rad. Assuming linear decline of the download factor until advance ratio
        0.5.

        Eq.:
            T sin(-alpha) = D + W sin(gamma)   (1) +
            T cos(-alpha) = W cos(gamma)       (2)
        """
        weight = self.mtow * self.gravity
        
        # Download factor
        if abs(advance_ratio) > 0.5: k_DL = 0
        else: k_DL = (1 - abs(advance_ratio) / 0.5) * self.download_factor

        # Thrust [N]
        return ((drag + weight * (np.sin(gamma) + np.cos(gamma))) 
                / (np.cos(alpha) - np.sin(alpha))) / (1 - k_DL)


    def get_angle_of_attack(self, drag, gamma):
        """ Calculate the angle of attack based on a force balance with thrust,
        weight, and fuselage drag. Climb angle gamma in rad.
        """
        weight = self.mtow * self.gravity

        # Angle of attack [rad]
        return - np.arctan(
            (drag + weight * np.sin(gamma)) 
            / (weight * np.cos(gamma)))


