# -*- coding: utf-8 -*-
"""
Created on 2022-01-30

@author: Fabian Grimm (f.grimm@tum.de)
"""

import numpy as np
import yaml
from modules.rotor import Rotor
from modules.mission import get_density

class Aircraft:
	"""
	"""
	def __init__(self, filename):
		# Load aircraft data from YAML file
		with open('configurations/' + filename) as file:
			data = yaml.safe_load(file)

		if 'Main rotor' in data:
			self.main_rotor = Rotor(data['Main rotor'])
		if 'Tail rotor' in data:
			self.tail_rotor = Rotor(data['Tail rotor'])

		self.power_available    = data['Engine']['Power available']
		self.sfc                = data['Engine']['SFC']
		self.download_factor    = data['Fuselage']['Download factor']
		self.drag_area          = data['Fuselage']['Drag area']
		self.empty_weight_ratio = data['Misc']['Empty weight ratio']
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
		self.density      = get_density(self.height, self.temp_offset)


	def get_parasite_power(self):
		""" Calculate the parasite power created by the fuselage drag in 
		forward flight.
		"""
		# Parasite power [W]
		return 0.5 * self.density * self.flight_speed ** 3 * self.drag_area


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


	def get_initial_mtow(self):
		""" Estimate the initial maximum take-off weight based on the mission 
		profile. [p.90-91]
		"""
		fuel_mass = self.sfc / 1e3 * self.power_available * self.duration
		useful_load = fuel_mass + self.payload + self.crew_mass

		# MTOW [kg]
		return useful_load / (1 - self.empty_weight_ratio)


	def iterate_first_sizing(self):
		""" Perform one step in the first sizing loop for conventional 
		helicopter configurations (requires main rotor).

		TODO: 
			- Solidity: Blade loading diagram instead of constant chord.
			- EW models [Prouty]
		"""
		# Main rotor design (more factors should be considered) [p.98-172]
		rotor = self.main_rotor
		rotor.solidity = rotor.get_solidity()
		weight = self.mtow * self.gravity
		rotor.radius = rotor.get_min_power_radius(self.density, thrust=weight)

		# Preliminary performance estimation [p.175]
		induced_power = rotor.get_induced_power(self.density, thrust=weight)
		profile_power = rotor.get_profile_power(self.density)
		self.hover_power = induced_power + profile_power

		# Mass estimation incl. fuel [p.177-179]
		empty_weight = self.empty_weight_ratio * self.mtow
		fuel_mass = self.sfc / 1e3 * 1.1 * self.hover_power * self.duration
		self.mtow = empty_weight + fuel_mass + self.payload + self.crew_mass

		# Disc loading and FM [p.181-183]
		weight = self.mtow * self.gravity
		rotor.disc_loading = rotor.get_disc_loading(thrust=weight)
		self.figure_of_merit = (induced_power / rotor.kappa) / self.hover_power


	def iterate_second_sizing(self):
		""" Perform one step in the second sizing loop for conventional 
		helicopter configurations (requires tail rotor).

		TODO: 
			- Solidity: Empirical reference instead of constant chord.
			- How should the download factor be considered? What counts as slow
				flight?
		"""
		# Tail rotor design [p.204-213]
		tr = self.tail_rotor
		tr.radius = 0.4 * np.sqrt(2.2 * self.mtow / 1000)
		tr.solidity = tr.get_solidity()

		# Fuselage[p.230]
		parasite_power = self.get_parasite_power()

		# Refined performance calculation
		drag = self.get_fuselage_drag()
		self.alpha = self.get_angle_of_attack(drag)
		thrust = self.get_thrust(drag)
		induced_velocity = self.main_rotor.get_induced_velocity(
			self.density, self.flight_speed, self.alpha, thrust)
		induced_power = self.main_rotor.kappa * thrust * induced_velocity

		# Selection of engine and gearbox
		# Mass estimation incl. fuel
	

	def preliminary_design(self, logs=True):
		""" Preliminary design.
		"""
		self.set_mission_segment(1)
		self.mtow = self.get_initial_mtow()

		if logs:
			print(f'Mission segment: 1')
			print(f'Height: {self.height:.2f} m')
			print(f'Density: {self.density:.2f} kg m^-3\n-----')
			print(f'First MTOW estimation: {self.mtow:10.2f} kg')

		self.iterate_first_sizing()
		self.iterate_second_sizing()

		if logs:
			print(f'Hover power: {(self.hover_power / 1000):20.2f} kW')
			print(f'Tail rotor radius: {self.tail_rotor.radius:14.2f} m')