# -*- coding: utf-8 -*-
"""
Created on 2022-02-08

@author: Fabian Grimm (f.grimm@tum.de)
"""

from aircraft import Aircraft
from preliminary_design import preliminary_design

"""
"""

# Create aircraft object with initial values from YAML file
aircraft = Aircraft('concept_01.yaml')

# Plot mission
aircraft.mission.plot_mission()

# Start premilinary design
preliminary_design(aircraft, logs=True)