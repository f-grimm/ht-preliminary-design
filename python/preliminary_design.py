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

# Add main and tail rotor
aircraft.main_rotor = Rotor('concept_01', 'Main rotor')
aircraft.tail_rotor = Rotor('concept_01', 'Tail rotor')

# Estimate MTOW based on mission profile
aircraft.mtow = aircraft.get_initial_mtow()
print(f'First MTOW estimation: {aircraft.mtow:7.1f} kg')

# Perform one step in the first sizing loop
aircraft.iterate_first_sizing()
print(f'Hover power: {(aircraft.hover_power / 1000):17.1f} kW')

# Perform one step in the second sizing loop
aircraft.iterate_second_sizing()



