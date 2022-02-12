# -*- coding: utf-8 -*-
"""
Created on 2022-02-08

@author: Fabian Grimm (f.grimm@tum.de)
"""

from modules.aircraft import Aircraft

"""
"""

# Create aircraft object with initial values from YAML file
aircraft = Aircraft('concept_01.yaml')

# Start premilinary design
aircraft.preliminary_design(logs=True)

