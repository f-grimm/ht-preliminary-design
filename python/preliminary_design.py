# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 19:36:45 2022

@author: Fabian Grimm (f.grimm@tum.de)
"""

import numpy as np

class Aircraft:
	"""
	"""
	def __init__(self):
		# TODO: Import as YAML
		# General
		self.mtow = 0 # kg
		self.main_rotor = Rotor()

		# Power
		self.power_available = 500 # kW

		# Mission
		self.payload   = 600 # kg
		self.crew_mass = 200 # kg
		self.duration  = 2.5 # h
		
		# Assumptions
		self.gravity         = 9.81 # ms^-2
		self.download_factor = 0.05 # -


	def get_initial_mtow(self):
		"""Estimate the initial maximum take-off mass.
		"""
		fuel_mass = 0.25 * self.power_available * self.duration
		useful_load = fuel_mass + self.payload + self.crew_mass
		empty_weight_ratio = 0.5

		# MTOW
		return useful_load / (1 - empty_weight_ratio)


	def iterate_first_sizing(self):
		"""Perform one step in the first sizing loop.
		"""
		# ---- Main rotor design ----------------------------------------------

		weight = self.mtow * self.gravity
		rotor = self.main_rotor
		rotor.solidity = rotor.get_solidity()
		rotor.radius = rotor.get_min_power_radius(thrust=weight)
		# (additional factors should be considered)

		# ---- Preliminary performance estimation -----------------------------
		
		induced_power = rotor.get_induced_power(thrust=weight)
		profile_power = rotor.get_profile_power()
		hover_power = induced_power + profile_power

		print('Hover power: {:17.1f} kW'.format(hover_power / 1000))

		# ---- Mass estimation incl. fuel -------------------------------------
		
		# Layton
		# self.mtow = ...


	def iterate_second_sizing(self):
		"""Perform one step in the second sizing loop.
		"""
		# Tail rotor design
		# Fuselage and accessories
		# Refined performance calculation
		# Selection of engine and gearbox
		# Mass estimation incl. fuel


class Rotor:
	"""
	"""
	def __init__(self):
		# Geometry
		self.radius           = 6.5 # m
		self.number_of_blades = 4 # -
		self.chord            = 0.27 # m
		self.solidity         = 1 # -
		
		# Aerodynamics
		self.density              = 1.225 # kgm^-3
		self.kappa                = 1.15 # -
		self.zero_lift_drag_coeff = 0.011 # -
		self.tip_velocity         = 213 # ms^-1


	def get_induced_power(self, thrust):
		"""Calculate the induced power.
		"""
		ideal_induced_power = np.sqrt(
			thrust ** 3 / (2 * self.density * np.pi * self.radius ** 2))
		
		# Induced power
		return self.kappa * ideal_induced_power


	def get_profile_power(self):
		"""Calculate the profile power.
		"""
		# Profile power
		return (1 / 8 * self.density * self.tip_velocity ** 3 * self.solidity 
				* self.zero_lift_drag_coeff * np.pi * self.radius ** 2)


	def get_min_power_radius(self, thrust):
		"""Calculate the optimal rotor radius w.r.t. induced and profile power 
		in hover.
		"""
		# Optimal radius
		return (1 / self.tip_velocity 
				* np.sqrt(2 * thrust / (self.density * np.pi)) 
		        * (self.kappa / (self.solidity * self.zero_lift_drag_coeff)) 
		          ** (1 / 3))


	def get_solidity(self):
		"""Calculate the rotor solidity (rectangular approximation).
		"""
		# Solidity
		return self.number_of_blades * self.chord / (np.pi * self.radius)	





if __name__ == "__main__":
	"""
	"""
	concept = Aircraft()

	concept.mtow = concept.get_initial_mtow()
	print("First MTOW estimation: {:7.1f} kg".format(concept.mtow))

	concept.iterate_first_sizing()

