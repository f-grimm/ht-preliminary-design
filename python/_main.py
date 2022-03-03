# -*- coding: utf-8 -*-
"""
Created on 2022-02-08

@author: Fabian Grimm (f.grimm@tum.de)
"""

from helicopter import Helicopter
from mission import Mission

"""
"""

# Create aircraft with initial values from YAML file
concept = Helicopter('concept_01')

# Load mission from YAML file
mission = Mission('hover_cruise') 

# Start premilinary design
concept.preliminary_design(mission, segment=1, logs=True)

# Plot mission profile
# mission.plot_mission()
