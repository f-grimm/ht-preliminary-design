# -*- coding: utf-8 -*-
"""
Created on 2022-02-08

@author: Fabian Grimm (f.grimm@tum.de)
"""

import numpy as np

class Mission():
	"""
	"""
	def __init__(self, mission_data: dict):
		
		self.data         = mission_data
		self.name         = mission_data['Name']
		self.payload      = self.data['Payload'][0]
		self.crew_mass    = self.data['Crew mass'][0]
		self.duration     = self.data['Duration'][0]
		self.flight_speed = self.data['Flight speed'][0]
		self.climb_angle  = self.data['Climb angle'][0] * np.pi / 180
		self.height       = self.data['Height'][0]
		self.temp_offset  = self.data['Temperature offset'][0]
		self.density      = self.get_density()


	def set_mission_segment(self, i):
		""" Select the mission segment from lists defined in the YAML file. 
		E.g. Payload: [600, 550] # kg
		       -> i =  0,   1
		"""
		self.payload      = self.data['Payload'][i]
		self.crew_mass    = self.data['Crew mass'][i]
		self.duration     = self.data['Duration'][i]
		self.flight_speed = self.data['Flight speed'][i]
		self.climb_angle  = self.data['Climb angle'][i] * np.pi / 180
		self.height       = self.data['Height'][i]
		self.temp_offset  = self.data['Temperature offset'][i]
		self.density      = self.get_density()

		
	def get_density(self):
		""" Calculate the density at a given height based on the international
		standard atmosphere (ISA). Deviation from the ISA is considered via a
		temperature offset. [p. 278]
		"""

		temperature_msl = 288.15
		pressure_msl = 101325
		radius_earth = 6356766
		gas_constant = 287.05     
		gamma = - 0.0065     
		n = 1.235
		
		geopotential_height = (radius_earth * self.height 
		                       / (radius_earth + self.height))
		temperature_isa = temperature_msl + gamma * geopotential_height
		temperature = temperature_isa + self.temp_offset
		pressure = (pressure_msl * (1 + gamma / temperature_msl * 
		            geopotential_height) ** (n / (n - 1)))

		# Density [kg m^-3]
		return pressure / (gas_constant * temperature)
