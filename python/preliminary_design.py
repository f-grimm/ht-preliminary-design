# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 19:36:45 2022

@author: Fabian Grimm (f.grimm@tum.de)
"""

import numpy as np

class Helicopter:
	"""
	"""
	def __init__(self, requirements, initial_parameters, assumptions):
		self.mtom      = 0
		self.solidity  = 1
		
		self.payload   = requirements["payload"]
		self.crew_mass = requirements["crew_mass"]
		self.duration  = requirements["duration"]

		self.radius           = initial_parameters["rotor_radius"]
		self.number_of_blades = initial_parameters["number_of_blades"]
		self.chord            = initial_parameters["chord"]
		self.power_available  = initial_parameters["power_available"]

		self.gravity              = assumptions["gravity"]
		self.density              = assumptions["density"]
		self.kappa                = assumptions["kappa"]
		self.zero_lift_drag_coeff = assumptions["zero_lift_drag_coeff"]
		self.tip_velocity         = assumptions["tip_velocity"]
		self.download_factor      = assumptions["download_factor"]


	def get_initial_mtom(self):
		"""Estimate the initial maximum take-off mass.
		"""
		P_av   = self.power_available
		t      = self.duration
		m_pl   = self.payload
		m_crew = self.crew_mass
		
		m_fuel = 0.25 * P_av * t
		useful_load = m_pl + m_fuel + m_crew
		empty_weight_ratio = 0.5

		# Return MTOM
		return useful_load / (1 - empty_weight_ratio)


	def get_induced_power(self, thrust):
		"""Calculate the induced power.
		"""
		T     = thrust
		rho   = self.density
		R     = self.radius
		kappa = self.kappa

		P_i_ideal = np.sqrt(T ** 3 / (2 * rho * np.pi * R ** 2))

		# Return induced power
		return kappa * P_i_ideal


	def get_profile_power(self):
		"""Calculate the profile power.
		"""
		rho    = self.density
		V_tip  = self.tip_velocity
		sigma  = self.solidity
		C_d0   = self.zero_lift_drag_coeff
		R      = self.radius

		# Return profile power
		return (1 / 8 * rho * V_tip ** 3 * sigma * C_d0 * np.pi * R ** 2)


	def get_min_power_radius(self, thrust):
		"""Calculate the optimal rotor radius w.r.t. induced and profile power 
		in hover.
		"""
		T     = thrust
		rho   = self.density
		V_tip = self.tip_velocity
		kappa = self.kappa
		sigma = self.solidity
		C_d0  = self.zero_lift_drag_coeff

		# Return optimal radius
		return (1 / V_tip * np.sqrt(2 * T / (rho * np.pi)) 
		        * (kappa / (sigma * C_d0)) ** (1 / 3))


	def get_solidity(self):
		"""Calculate the rotor solidity (rectangular approximation).
		"""
		N_b = self.number_of_blades
		c   = self.chord
		R   = self.radius

		# Return solidity
		return N_b * c / (np.pi * R)


	def iterate_first_sizing(self):
		"""Perform one step in the first sizing loop.
		"""
		# ---- Main rotor design ----------------------------------------------

		weight = self.mtom * self.gravity
		self.solidity = self.get_solidity()
		self.radius = self.get_min_power_radius(thrust=weight)
		# (additional factors should be considered)

		# ---- Preliminary performance estimation -----------------------------
		
		P_i = self.get_induced_power(thrust=weight)
		P_0 = self.get_profile_power()
		P_h = P_i + P_0

		# ---- Mass estimation incl. fuel -------------------------------------
		
		# Layton
		# self.mtom = ...


	def iterate_second_sizing(self):
		"""Perform one step in the second sizing loop.
		"""
		# Tail rotor design
		# Fuselage and accessories
		# Refined performance calculation
		# Selection of engine and gearbox
		# Mass estimation incl. fuel


	def loop(self):
		"""
		"""
		converged = False

		while not converged:
			self.iterate_first_sizing()
			
			if converged:
				self.iterate_second_sizing()


		
if __name__ == "__main__":
	"""
	"""
	requirements = {
		"payload":   600, # kg
		"crew_mass": 200, # kg
		"duration":  2.5, # h
		}

	initial_parameters = {
		"rotor_radius":     6.5, # m
		"number_of_blades": 4, # -
		"chord":            0.27, # m
		"power_available":  500, # kW
		}

	assumptions = {
		"gravity":              9.81, # ms^-2
		"density":              1.225, # kgm^-3
		"kappa":                1.15, # -
		"zero_lift_drag_coeff": 0.011, # -
		"tip_velocity":         213, # ms^-1
		"download_factor":      0.05, # -
		}

	concept = Helicopter(requirements, initial_parameters, assumptions)

	concept.mtom = concept.get_initial_mtom()
	print("First MTOM estimation: {} kg \n".format(concept.mtom))

	
	concept.iterate_first_sizing()

