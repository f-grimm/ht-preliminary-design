# -*- coding: utf-8 -*-
"""
Created on 2022-01-30

@author: Fabian Grimm (f.grimm@tum.de)
"""

import numpy as np

class Rotor:
	"""
	"""
	def __init__(self, rotor_data: dict):
		
		self.radius               = rotor_data['Radius']
		self.number_of_blades     = rotor_data['Number of blades']
		self.chord                = rotor_data['Chord']
		self.solidity             = rotor_data['Solidity']
		self.kappa                = rotor_data['Kappa']
		self.zero_lift_drag_coeff = rotor_data['Zero-lift drag coeff.']
		self.tip_velocity         = rotor_data['Tip velocity']

		if 'Power fraction' in rotor_data:
			self.power_fraction = rotor_data['Power fraction']


	def get_solidity(self):
		""" Calculate the rotor solidity (rectangular approximation).
		"""
		# Solidity [-]
		return self.number_of_blades * self.chord / (np.pi * self.radius)


	def get_disc_loading(self, thrust):
		""" Calculate the disc loading.
		"""
		# Disc loading [N m^-2]
		return thrust / (np.pi * self.radius ** 2)


	def get_induced_power_hover(self, density, thrust):
		""" Calculate the induced power in hover.
		"""
		ideal_induced_power = np.sqrt(
			thrust ** 3 / (2 * density * np.pi * self.radius ** 2))
		
		# Induced power [W]
		return self.kappa * ideal_induced_power


	def get_induced_velocity(self, density, v_inf, alpha, thrust):
		""" Calculate the induced velocity iteratively in forward flight (not 
		valid for low sink rate in axial flight) [p. 253]
		"""
		area = np.pi * self.radius ** 2
		v_i = np.sqrt(thrust / (2 * density * area))
		counter = 0
		error = 1

		while error > 1e-4:
			combined_velocity = np.sqrt(
				v_inf ** 2 - 2 * v_inf * v_i * np.sin(alpha) + v_i ** 2)
			v_i_new = thrust / (2 * density * area * combined_velocity)
			error = abs((v_i_new - v_i) / v_i)
			v_i = v_i_new
			counter += 1
			if counter > 1e6: raise ValueError('No convergence.')

		# Induced velocity [m s^-1]
		return v_i


	def get_induced_velocity_level(self, density, v_inf, thrust):
		""" Calculate the induced velocity using the approximate solution for 
		level flight (W >> D). [p. 254]
		"""
		v_h = np.sqrt(thrust / (2 * density * np.pi * self.radius ** 2))

		# Induced velocity [m s^-1]
		return np.sqrt(-0.5 * v_inf ** 2 
		               + np.sqrt(v_h ** 4 + 0.25 * v_inf ** 4))


	def get_profile_power(self, density, advance_ratio):
		""" Calculate the profile power. [p.258]
		"""
		# Profile power [W]
		return (1 / 8 * density * self.tip_velocity ** 3 * self.solidity 
		        * self.zero_lift_drag_coeff * np.pi * self.radius ** 2
		        * (1 + 4.65 * advance_ratio ** 2))


	def get_min_power_radius(self, density, thrust):
		""" Calculate the optimal rotor radius w.r.t. induced and profile power 
		in hover. [p.103]
		"""
		# Optimal radius [m]
		return (1 / self.tip_velocity 
		        * np.sqrt(2 * thrust / (density * np.pi)) 
		        * (self.kappa / (self.solidity * self.zero_lift_drag_coeff)) 
		          ** (1 / 3))


	def in_ground_effect(self, induced_power_oge):
		""" Calculate the required hover power in ground effect (IGE) according 
		to Hayden. [p. 288]
		"""
		k_G = 1 / (0.9926 + 0.0379 * (2 * self.radius / rotor_height) ** 2)

		# Power in ground effect [W]
		return k_G * induced_power_oge
