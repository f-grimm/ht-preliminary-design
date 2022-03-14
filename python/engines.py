# -*- coding: utf-8 -*-
"""
Created on 2022-02-08

@author: Fabian Grimm (f.grimm@tum.de)
"""

import numpy as np

class Engines():
	"""
	"""
	def __init__(self, engine_data:dict):
		# Get values from engine data dict (default None if key is unavailable)
		self.number_of_engines = engine_data.get('Number of engines', 1)
		self.power_available   = engine_data.get('Power available')
		self.sfc               = engine_data.get('SFC')
		self.a                 = engine_data.get('A')
		self.b                 = engine_data.get('B')


	def get_sfc(self, temperature_ratio, pressure_ratio, power):
		""" Calculate the specific fuel consumption based on the engine 
		parameters A and B. [p.310]
		"""
		# SFC [kg/Wh]
		return (self.a * pressure_ratio * np.sqrt(temperature_ratio) / power
				+ self.b)


