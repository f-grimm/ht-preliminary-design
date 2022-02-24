# -*- coding: utf-8 -*-
"""
Created on 2022-02-08

@author: Fabian Grimm (f.grimm@tum.de)
"""

import numpy as np
import matplotlib.pyplot as plt

class Mission():
	"""
	"""
	def __init__(self, mission_data: dict):

		self.name = mission_data['Name']
		self.data = mission_data['Mission']
		self.set_mission_segment(0)
		

	def set_mission_segment(self, i):
		""" Select the mission segment from lists defined in the YAML file. 
		E.g. Payload: [600, 550] # kg
		       -> i =  0,   1
		"""
		self.segment      = i
		self.duration     = self.data['Duration'][i]
		self.payload      = self.data['Payload'][i]
		self.crew_mass    = self.data['Crew mass'][i]
		self.climb_angle  = self.data['Climb angle'][i] * np.pi / 180
		self.flight_speed = self.data['Flight speed'][i]
		self.temp_offset  = self.data['Temperature offset'][i]
		self.height       = self.data['Height'][i]
		(self.density, 
		self.pressure, 
		self.temperature) = self.atmosphere()

		
	def atmosphere(self):
		""" Calculate the density at a given height based on the international
		standard atmosphere (ISA). Deviation from the ISA is considered via a
		temperature offset. [p. 278]
		"""
		temperature_msl = 288.15
		pressure_msl = 101325
		radius_earth = 6356766
		gas_constant = 287.05     
		gamma = - 0.0065     
		n = 1.235
		
		geopotential_height = (radius_earth * self.height 
		                       / (radius_earth + self.height))
		temperature_isa = temperature_msl + gamma * geopotential_height
		temperature = temperature_isa + self.temp_offset
		pressure = (pressure_msl * (1 + gamma / temperature_msl * 
		            geopotential_height) ** (n / (n - 1)))
		density = pressure / (gas_constant * temperature)

		# Density [kg m^-3], pressure [Pa], temperature [K]
		return density, pressure, temperature


	def plot_mission(self):
		""" Plot height, payload, and crew mass over the duration of the 
		mission.
		"""
		# Figure settings
		fig, ax = plt.subplots(figsize=(8,5))
		ax.set_title('Mission: ' + self.name, fontweight='bold')
		ax_1 = ax.twinx()
		ax.grid(True)
		ax_1.grid(False)

		# Time list with double elements between segments to plot constant 
		# parameters, e.g. [0, 2.5, 2.5, 3.5]
		time = ([0] + [t for t in np.cumsum(self.data['Duration']) 
		               for _ in (0, 1)][:-1])

		# Constant parameters
		payload = [m for m in self.data['Payload'] for _ in (0, 1)]
		
		# Linear parameters
		height = [h for h in self.data['Height'] for _ in (0, 1)][1:-1]

		# Plot, axes
		ax.set_xlabel('Time [h]')
		ax.set_ylabel('Height [m]')
		p, = ax.plot(time, height, color='tab:red', label='Height')
		ax.set_ylim(bottom=0)
		ax_1.set_ylabel('Mass [kg]')
		p_1, = ax_1.plot(time, payload, color='tab:blue', label='Payload')
		ax_1.set_ylim(bottom=0, top=max(payload) / 0.5)
		
		# Legend
		plots = [p, p_1]
		legend = ax.legend(handles=plots, loc='best')

		# Layout
		fig.tight_layout()
		plt.show()
