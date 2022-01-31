# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 19:36:45 2022

@author: Fabian Grimm (f.grimm@tum.de)
"""

from modules.aircraft import *
from modules.rotor import *

"""
"""

# Create aircraft object with initial values from YAML file
aircraft = Aircraft('concept_01')

# Add main rotor
aircraft.main_rotor = Rotor('concept_01')

# Estimate MTOW based on mission profile
aircraft.mtow = aircraft.get_initial_mtow()
print("First MTOW estimation: {:7.1f} kg".format(aircraft.mtow))

# Perform one step in the first sizing loop
aircraft.iterate_first_sizing()
print('Hover power: {:17.1f} kW'.format(aircraft.hover_power / 1000))


