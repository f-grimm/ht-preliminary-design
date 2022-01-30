# -*- coding: utf-8 -*-
"""
Created on 2022-01-30

@author: Fabian Grimm (f.grimm@tum.de)
"""

import numpy as np

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
		
		# Induced power [W]
		return self.kappa * ideal_induced_power


	def get_profile_power(self):
		"""Calculate the profile power.
		"""
		# Profile power [W]
		return (1 / 8 * self.density * self.tip_velocity ** 3 * self.solidity 
				* self.zero_lift_drag_coeff * np.pi * self.radius ** 2)


	def get_min_power_radius(self, thrust):
		"""Calculate the optimal rotor radius w.r.t. induced and profile power 
		in hover.
		"""
		# Optimal radius [m]
		return (1 / self.tip_velocity 
				* np.sqrt(2 * thrust / (self.density * np.pi)) 
		        * (self.kappa / (self.solidity * self.zero_lift_drag_coeff)) 
		          ** (1 / 3))


	def get_solidity(self):
		"""Calculate the rotor solidity (rectangular approximation).
		"""
		# Solidity [-]
		return self.number_of_blades * self.chord / (np.pi * self.radius)