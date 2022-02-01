# -*- coding: utf-8 -*-
"""
Created on 2022-01-30

@author: Fabian Grimm (f.grimm@tum.de)
"""

import numpy as np
import yaml

class Aircraft:
	"""
	"""
	def __init__(self, filename):
		# Load aircraft data from YAML file
		with open('configurations/' + filename + '.yaml') as file:
			aircraft_data = yaml.safe_load(file)['Aircraft']

		self.mtow            = aircraft_data['MTOW']
		self.power_available = aircraft_data['Power available']
		self.hover_power     = aircraft_data['Hover power']
		self.payload         = aircraft_data['Payload']
		self.crew_mass       = aircraft_data['Crew mass']
		self.duration        = aircraft_data['Duration']
		self.flight_speed    = aircraft_data['Flight speed']
		self.gravity         = aircraft_data['Gravity']
		self.download_factor = aircraft_data['Download factor']
		self.drag_area       = aircraft_data['Drag area']


	def get_initial_mtow(self):
		""" Estimate the initial maximum take-off weight based on the mission 
		profile. [p.90-91]
		"""
		fuel_mass = 0.25e-3 * self.power_available * self.duration
		useful_load = fuel_mass + self.payload + self.crew_mass
		empty_weight_ratio = 0.5

		# MTOW [kg]
		return useful_load / (1 - empty_weight_ratio)


	def get_parasite_power(self, density):
		""" Calculate the parasite power created by the fuselage drag in 
		forward flight.
		"""
		# Parasite power [W]
		return 0.5 * density * self.flight_speed ** 3 * self.drag_area


	def iterate_first_sizing_conv(self):
		""" Perform one step in the first sizing loop for conventional 
		helicopter configurations (Requires self.main_rotor of the class 
		"Rotor").

		TODO: Solidity: Blade loading diagram instead of constant chord.
		"""
		# Main rotor design (more factors should be considered) [p.98-172]
		rotor = self.main_rotor
		rotor.solidity = rotor.get_solidity()
		weight = self.mtow * self.gravity
		rotor.radius = rotor.get_min_power_radius(thrust=weight)

		# Preliminary performance estimation [p.175]
		induced_power = rotor.get_induced_power(thrust=weight)
		profile_power = rotor.get_profile_power()
		self.hover_power = induced_power + profile_power

		# Mass estimation incl. fuel [p.177-179]
		empty_weight = 1000 # (TODO: Find empirical models)
		fuel_mass = 0.38e-3 * 1.1 * self.hover_power * self.duration
		self.mtow = empty_weight + fuel_mass + self.payload + self.crew_mass

		# Disc loading and FM [p.181-183]
		weight = self.mtow * self.gravity
		rotor.disc_loading = rotor.get_disc_loading(thrust=weight)
		self.figure_of_merit = (induced_power / rotor.kappa) / self.hover_power


	def iterate_second_sizing_conv(self):
		""" Perform one step in the second sizing loop.

		TODO: 
			- Solidity: Empirical reference instead of constant chord.
			- How should the download factor be considered? What counts as slow
				flight?
			- Altitude and ISA calculation instead of density in rotor class
		"""
		# Tail rotor design [p.204-213]
		tr = self.tail_rotor
		tr.radius = 0.4 * np.sqrt(2.2 * self.mtow / 1000)
		tr.solidity = tr.get_solidity()

		# Fuselage and accessories [p.228]
		parasite_power = get_parasite_power(tr.density)


		# Refined performance calculation
		# Selection of engine and gearbox
		# Mass estimation incl. fuel