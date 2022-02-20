# -*- coding: utf-8 -*-
"""
Created on 2022-02-08

@author: Fabian Grimm (f.grimm@tum.de)
"""

from aircraft import Aircraft
from preliminary_design import preliminary_design, plot_powers

"""
"""

# Create aircraft with initial values from YAML file
aircraft = Aircraft('concept_01', mission='hc')

# Start premilinary design
preliminary_design(aircraft, mission_segment=1, logs=True)

# Plot mission profile
# aircraft.mission.plot_mission()

# Plot power over flight speed
# plot_powers(aircraft, max_velocity=80)
