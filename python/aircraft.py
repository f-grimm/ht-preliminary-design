# -*- coding: utf-8 -*-
"""
Created on 2022-01-30

@author: Fabian Grimm (f.grimm@tum.de)
"""

import numpy as np
import yaml
import rotor
import mission

class Aircraft:
	"""
	"""
	def __init__(self, filename):
		# Load aircraft data from YAML file
		with open('configurations/' + filename) as file:
			data = yaml.safe_load(file)

		if 'Main rotor' in data:
			self.main_rotor = rotor.Rotor(data['Main rotor'])
		if 'Tail rotor' in data:
			self.tail_rotor = rotor.Rotor(data['Tail rotor'])

		self.power_available    = data['Engine']['Power available']
		self.accessory_power    = data['Engine']['Accessory power']
		self.sfc                = data['Engine']['SFC'] / 1000
		self.download_factor    = data['Fuselage']['Download factor']
		self.drag_area          = data['Fuselage']['Drag area']
		self.empty_weight_ratio = data['Misc']['Empty weight ratio']
		self.eta_transmission   = data['Misc']['Transmission efficiency']
		self.gravity            = data['Misc']['Gravity']
		self.mission            = data['Mission']


	def set_mission_segment(self, i):
		""" Select the mission segment from lists defined in the YAML file. 
		E.g. Payload: [600, 550] # kg
		       -> i =  1,   2
		"""
		self.payload      = self.mission['Payload'][i - 1]
		self.crew_mass    = self.mission['Crew mass'][i - 1]
		self.duration     = self.mission['Duration'][i - 1]
		self.flight_speed = self.mission['Flight speed'][i - 1]
		self.climb_angle  = self.mission['Climb angle'][i - 1] * np.pi / 180
		self.height       = self.mission['Height'][i - 1]
		self.temp_offset  = self.mission['Temperature offset'][i - 1]
		self.density      = mission.get_density(self.height, self.temp_offset)


	def get_parasite_power(self):
		""" Calculate the parasite power created by the fuselage drag in 
		forward flight.
		"""
		# Parasite power [W]
		return 0.5 * self.density * self.flight_speed ** 3 * self.drag_area


	def get_climb_power(self):
		""" Claculate the climb power due to change in potential energy.
		"""
		weight = self.mtow * self.gravity

		# Climb power [W]
		return weight * self.flight_speed * np.sin(self.climb_angle)


	def get_fuselage_drag(self):
		""" Calculate the fuselage drag.
		"""
		# Fuselage drag [N]
		return 0.5 * self.density * self.flight_speed ** 2 * self.drag_area


	def get_thrust(self, drag):
		""" Calculate the required thrust in translatory motion.
		"""
		weight = self.mtow * self.gravity

		# Thrust [N]
		return ((drag + weight * (np.sin(self.climb_angle) 
		         + np.cos(self.climb_angle))) / 
		        (np.sin(self.alpha) + np.cos(self.alpha)))


	def get_angle_of_attack(self, drag):
		""" Calculate the angle of attack based on a force balance with thrust,
		weight, and fuselage drag.
		"""
		weight = self.mtow * self.gravity

		# Angle of attack [rad]
		return np.arctan((drag + weight * np.sin(self.climb_angle)) /
		                 (weight * np.cos(self.climb_angle)))