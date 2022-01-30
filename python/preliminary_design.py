# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 19:36:45 2022

@author: Fabian Grimm (f.grimm@tum.de)
"""

from modules.aircraft import *
from modules.rotor import *


concept = Aircraft()
concept.main_rotor = Rotor()

concept.mtow = concept.get_initial_mtow()
print("First MTOW estimation: {:7.1f} kg".format(concept.mtow))

concept.iterate_first_sizing()
print('Hover power: {:17.1f} kW'.format(concept.hover_power / 1000))


