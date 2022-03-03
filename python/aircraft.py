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
		self.power_available    = data['Engine']['Power available']
		self.accessory_power    = data['Engine']['Accessory power']
		self.sfc                = data['Engine']['SFC'] * 1e-3
		self.download_factor    = data['Fuselage']['Download factor']
		self.drag_area          = data['Fuselage']['Drag area']
		self.landing_gear_type  = data['Landing gear']['Type']
		self.number_of_legs     = data['Landing gear']['Number of legs']
		self.empty_weight_ratio = data['Misc']['Empty weight ratio']
		self.eta_transmission   = data['Misc']['Transmission efficiency']
		self.special_equipment  = data['Misc']['Special equipment']
		self.number_of_seats    = data['Misc']['Number of seats']
		self.gravity            = data['Misc']['Gravity']
		self.mtow               = 0
		self.alpha              = 0
		self.drag               = 0


	def get_parasite_power(self, density, flight_speed):
		""" Calculate the parasite power created by the fuselage drag in 
		forward flight.
		"""
		# Parasite power [W]
		return (0.5 * density * flight_speed ** 3 * self.drag_area)


	def get_climb_power(self, flight_speed, climb_angle):
		""" Claculate the climb power due to change in potential energy.
		"""
		weight = self.mtow * self.gravity

		# Climb power [W]
		return (weight * flight_speed * np.sin(climb_angle))


	def get_fuselage_drag(self, density, flight_speed):
		""" Calculate the fuselage drag.
		"""
		# Fuselage drag [N]
		return (0.5 * density * flight_speed ** 2 * self.drag_area)


	def get_thrust(self, climb_angle):
		""" Calculate the required thrust in translatory motion.
		Eq.:
			T sin(-alpha) = D + W sin(gamma)   (1) +
			T cos(-alpha) = W cos(gamma)       (2)
		"""
		weight = self.mtow * self.gravity

		# Thrust [N]
		return ((self.drag + weight * (np.sin(climb_angle) 
		        + np.cos(climb_angle))) 
		        / (np.cos(self.alpha) - np.sin(self.alpha)))


	def get_angle_of_attack(self, climb_angle):
		""" Calculate the angle of attack based on a force balance with thrust,
		weight, and fuselage drag.
		"""
		weight = self.mtow * self.gravity

		# Angle of attack [rad]
		return - np.arctan(
			(self.drag + weight * np.sin(climb_angle)) 
			/ (weight * np.cos(climb_angle)))


