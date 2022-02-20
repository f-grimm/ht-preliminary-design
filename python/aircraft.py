# -*- coding: utf-8 -*-
"""
Created on 2022-01-30

@author: Fabian Grimm (f.grimm@tum.de)
"""

import numpy as np
import yaml
import rotor
import mission as mission_module

class Aircraft:
	"""
	"""
	def __init__(self, filename: str, mission: str):
		# Load aircraft data from YAML file
		with open('configurations/' + filename + '.yaml') as file:
			data = yaml.safe_load(file)

		if 'Main rotor' in data:
			self.main_rotor = rotor.Rotor(data['Main rotor'])
		if 'Tail rotor' in data:
			self.tail_rotor = rotor.Rotor(data['Tail rotor'])

		self.name               = data['Name']
		self.power_available    = data['Engine']['Power available']
		self.accessory_power    = data['Engine']['Accessory power']
		self.sfc                = data['Engine']['SFC'] / 1000
		self.download_factor    = data['Fuselage']['Download factor']
		self.drag_area          = data['Fuselage']['Drag area']
		self.empty_weight_ratio = data['Misc']['Empty weight ratio']
		self.eta_transmission   = data['Misc']['Transmission efficiency']
		self.gravity            = data['Misc']['Gravity']
		self.mtow               = 0
		self.alpha              = 0
		self.drag               = 0

		# Load mission data from YAML file
		with open('configurations/missions/' + mission + '.yaml') as file:
			mission_data = yaml.safe_load(file)

		self.mission = mission_module.Mission(mission_data)


	def get_parasite_power(self):
		""" Calculate the parasite power created by the fuselage drag in 
		forward flight.
		"""
		# Parasite power [W]
		return (0.5 * self.mission.density * self.mission.flight_speed ** 3 
		        * self.drag_area)


	def get_climb_power(self):
		""" Claculate the climb power due to change in potential energy.
		"""
		weight = self.mtow * self.gravity

		# Climb power [W]
		return (weight * self.mission.flight_speed 
		        * np.sin(self.mission.climb_angle))


	def get_fuselage_drag(self):
		""" Calculate the fuselage drag.
		"""
		# Fuselage drag [N]
		return (0.5 * self.mission.density * self.mission.flight_speed ** 2 
		        * self.drag_area)


	def get_thrust(self):
		""" Calculate the required thrust in translatory motion.
		Eq.:
			T sin(-alpha) = D + W sin(gamma)   (1) +
			T cos(-alpha) = W cos(gamma)       (2)
		"""
		weight = self.mtow * self.gravity

		# Thrust [N]
		return ((self.drag + weight * (np.sin(self.mission.climb_angle) 
		         + np.cos(self.mission.climb_angle))) 
		        / (np.cos(self.alpha) - np.sin(self.alpha)))


	def get_angle_of_attack(self):
		""" Calculate the angle of attack based on a force balance with thrust,
		weight, and fuselage drag.
		"""
		weight = self.mtow * self.gravity

		# Angle of attack [rad]
		return - np.arctan(
			(self.drag + weight * np.sin(self.mission.climb_angle)) 
			/ (weight * np.cos(self.mission.climb_angle)))