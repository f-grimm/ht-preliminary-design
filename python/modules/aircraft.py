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
		self.gravity         = aircraft_data['Gravity']
		self.download_factor = aircraft_data['Download factor']


	def get_initial_mtow(self):
		""" Estimate the initial maximum take-off weight.
		"""
		fuel_mass = 0.25e-3 * self.power_available * self.duration
		useful_load = fuel_mass + self.payload + self.crew_mass
		empty_weight_ratio = 0.5

		# MTOW [kg]
		return useful_load / (1 - empty_weight_ratio)


	def iterate_first_sizing(self):
		""" Perform one step in the first sizing loop. 
		(Requires self.main_rotor of the class "Rotor")
		"""
		# Main rotor design (additional factors should be considered)
		weight = self.mtow * self.gravity
		rotor = self.main_rotor
		rotor.solidity = rotor.get_solidity()
		rotor.radius = rotor.get_min_power_radius(thrust=weight)

		# Preliminary performance estimation
		induced_power = rotor.get_induced_power(thrust=weight)
		profile_power = rotor.get_profile_power()
		self.hover_power = induced_power + profile_power

		# Mass estimation incl. fuel
		# Layton
		# self.mtow = ...


	def iterate_second_sizing(self):
		""" Perform one step in the second sizing loop.
		"""
		# Tail rotor design
		# Fuselage and accessories
		# Refined performance calculation
		# Selection of engine and gearbox
		# Mass estimation incl. fuel