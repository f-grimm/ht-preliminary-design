# -*- coding: utf-8 -*-
"""
Created on 2022-02-08

@author: Fabian Grimm (f.grimm@tum.de)
"""

def get_density(height, temperature_offset):
	""" Calculate the density at a given height based on the international
	standard atmosphere (ISA). Deviation from the ISA is considered via a
	temperature offset.
	"""

	temperature_msl = 288.15
	pressure_msl = 101325
	radius_earth = 6356766
	gas_constant = 287.05     
	gamma = - 0.0065     
	n = 1.235
	
	geopotential_height = radius_earth * height / (radius_earth + height)
	temperature_isa = temperature_msl + gamma * geopotential_height
	temperature = temperature_isa + temperature_offset
	pressure = (pressure_msl * (1 + gamma / temperature_msl * 
	            geopotential_height) ** (n / (n - 1)))

	# Density [kg m^-3]
	return pressure / (gas_constant * temperature)