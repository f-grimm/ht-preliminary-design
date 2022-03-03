# -*- coding: utf-8 -*-
"""
Created on 2022-02-12

@author: Fabian Grimm (f.grimm@tum.de)
"""

from aircraft import Aircraft
from mission import Mission
from rotor import Rotor
import yaml
import numpy as np
import matplotlib.pyplot as plt

class Helicopter(Aircraft):
	""" Aircraft sub-class for conventional helicopter configurations.

	Requires:
		Main rotor
		Tail rotor
	"""
	def __init__(self, filename: str):		
		Aircraft.__init__(self, filename)

		# Load data from YAML file
		with open('data/configurations/' + filename + '.yaml') as file:
			data = yaml.safe_load(file)

		self.main_rotor = Rotor(data['Main rotor'])
		self.tail_rotor = Rotor(data['Tail rotor'])


	def preliminary_design(self, mission: Mission, segment: int, logs=True):
		"""
		"""
		mission.set_mission_segment(segment)

		if logs:
			print('\nPreliminary Design\n' + '-' * 18)
			print(f'\nAicraft: {self.name}')
			print(f'Mission: {mission.name}')
			print(f'Segment: {mission.segment} '
			      + f'({mission.height:.0f}m, '
			      + f'{mission.duration:.2f}h, '
			      + f'{mission.flight_speed:.0f}m/s)')
			print(f'\nEstimate initial MTOW'
			      + f'\n(SFC {self.sfc * 1e3:.2f}, '
			      + f'EWR {self.empty_weight_ratio:.2f})')

		self.mtow = self.get_initial_mtow(mission)
			
		if logs:
			print(f'\n - MTOW: {self.mtow:27.2f} kg')
			print('\nEnter design loop')

		mtow_list = self.separate_loops(mission)

		if logs:
			print(f'Iterations: {len(mtow_list) - 1}')
			print(f'\n - Main rotor radius: {self.main_rotor.radius:14.2f} m')
			print(f' - Tail rotor radius: {self.tail_rotor.radius:14.2f} m')
			print(f' - Hover power: {(self.hover_power * 1e-3):20.2f} kW')
			print(f' - Drag: {self.drag:27.2f} N')
			print(f' - Thrust: {self.thrust:25.2f} N')
			print(f' - Angle of attack: {self.alpha * 180 / np.pi:16.2f} deg')
			print(f' - Advance ratio: {self.advance_ratio:18.2f} -')
			print(' - Induced velocity: '
			      + f'{self.main_rotor.induced_velocity:15.2f} m/s')
			print(f' - Total power: {(self.total_power * 1e-3):20.2f} kW')
			print(f' - MTOW: {self.mtow:27.2f} kg')
			print('')

			self.plot_mtow_convergence(mtow_list)


	def get_initial_mtow(self, mission: Mission):
		""" Estimate the initial maximum take-off weight based on the mission 
		profile; (SFC and empty weight ratio are constant). [pp.90-91]
		"""
		fuel_mass = (self.sfc * self.power_available * mission.duration)
		useful_load = (fuel_mass + mission.payload + mission.crew_mass)

		# MTOW [kg]
		return useful_load / (1 - self.empty_weight_ratio)


	def design_loop(self, mission: Mission):
		""" Combine first sizing as an inner loop with the second sizing in the 
		outer loop in order to determine the maximum take-off weight. [p.86]
		"""
		mtow_list = [self.mtow]
		counter = 0
		error = 1

		while error > 1e-4:

			while error > 0.1:
				self.iterate_first_sizing(mission)
				mtow_list.append(self.mtow)
				error = abs((mtow_list[-1] - mtow_list[-2]) / mtow_list[-1])
				counter += 1
				if counter > 1e3: break

			self.iterate_second_sizing(mission)
			mtow_list.append(self.mtow)
			error = abs((mtow_list[-1] - mtow_list[-2]) / mtow_list[-1])
			counter += 1
			if counter > 1e3: break

		# List of MTOW values [kg]
		return mtow_list


	def separate_loops(self, mission: Mission):
		""" Combine first sizing and second sizing as separate loops in order to 
		determine the maximum take-off weight.
		"""
		mtow_list = [self.mtow]
		counter = 0
		error = 1

		while error > 0.1:
			self.iterate_first_sizing(mission)
			mtow_list.append(self.mtow)
			error = abs((mtow_list[-1] - mtow_list[-2]) / mtow_list[-1])
			counter += 1
			if counter > 1e3: break

		counter = 0
		error = 1

		while error > 1e-4:
			self.iterate_second_sizing(mission)
			mtow_list.append(self.mtow)
			error = abs((mtow_list[-1] - mtow_list[-2]) / mtow_list[-1])
			counter += 1
			if counter > 1e3: break

		# List of MTOW values [kg]
		return mtow_list


	def iterate_first_sizing(self, mission):
		""" Perform one step in the first sizing loop [p.86]

		TODO: 
			Solidity: Blade loading diagram instead of constant chord.
		"""
		# ---- Main rotor sizing [pp.98-172] --------------------------------------

		# Rotor radius for min. power in hover (more factors should be considered)
		self.main_rotor.solidity = self.main_rotor.get_solidity()
		weight = self.mtow * self.gravity
		self.main_rotor.radius = self.main_rotor.get_min_power_radius(
			mission.density, thrust=weight)

		# ---- Preliminary performance estimation [p.175] -------------------------

		# Hover power
		induced_power = self.main_rotor.get_induced_power_hover(
			mission.density, thrust=weight)
		profile_power = self.main_rotor.get_profile_power(
			mission.density, advance_ratio=0)
		self.hover_power = induced_power + profile_power

		# ---- Mass estimation incl. fuel [pp.177-179] ----------------------------

		# Fuel consumption in hover, empty weight
		self.fuel_mass = (self.sfc * 1.1 * self.hover_power 
		                      * mission.duration)
		self.empty_weight = sum(self.mass_estimation(
			power=self.hover_power).values())
		self.mtow = (self.empty_weight + self.fuel_mass 
		                 + mission.payload + mission.crew_mass)

		# Disc loading and FM [pp.181-183]
		weight = self.mtow * self.gravity
		self.main_rotor.disc_loading = self.main_rotor.get_disc_loading(
			thrust=weight)
		self.main_rotor.figure_of_merit = (
			(induced_power / self.main_rotor.kappa) / self.hover_power)


	def iterate_second_sizing(self, mission):
		""" Perform one step in the second sizing loop. [p.86]

		TODO: 
			Solidity: Empirical reference instead of constant chord.
			How should the download factor be considered? What counts as slow
				flight?
			Engine calc., SFC = f(P) [p.310]
			Conditions at the beginning of each segment assumed. Half altitude for 
				climb?
		"""
		# ---- Tail rotor design [pp.204-213] -------------------------------------

		# Tail rotor design according to Layton
		self.tail_rotor.radius = 0.4 * np.sqrt(2.2 * self.mtow * 1e-3)
		self.tail_rotor.solidity = self.tail_rotor.get_solidity()

		# ---- Refined performance calculation [pp.251-321] -----------------------

		# Flight state
		self.drag = self.get_fuselage_drag(mission.density, mission.flight_speed)
		self.alpha = self.get_angle_of_attack(mission.climb_angle)
		self.thrust = self.get_thrust(mission.climb_angle)
		self.advance_ratio = (mission.flight_speed 
		                          / self.main_rotor.tip_velocity)
		self.main_rotor.induced_velocity = (
			self.main_rotor.get_induced_velocity(
				mission.density, mission.flight_speed, 
				self.alpha, self.thrust))

		# Power calculation
		induced_power = (self.main_rotor.kappa * self.thrust 
		                 * self.main_rotor.induced_velocity)
		profile_power = self.main_rotor.get_profile_power(
			mission.density, self.advance_ratio)
		parasite_power = self.get_parasite_power(mission.density, mission.flight_speed)
		climb_power = self.get_climb_power(mission.flight_speed, mission.climb_angle)
		main_rotor_power = (induced_power + profile_power + parasite_power 
		                    + climb_power)
		tail_rotor_power = self.tail_rotor.power_fraction * main_rotor_power
		transmission_losses = (((1 / self.eta_transmission) - 1) 
		                       * main_rotor_power)
		self.total_power = (main_rotor_power + tail_rotor_power 
		                        + transmission_losses + self.accessory_power)

		# Engine requirement
		temperature_ratio = mission.temperature / 288.15
		pressure_ratio = mission.pressure / 101325
		power_msl = (self.total_power * np.sqrt(temperature_ratio) 
		             / pressure_ratio)

		# ---- Refined mass estimation [pp.324-325] -------------------------------
		
		# Fuel mass, empty weight
		self.fuel_mass = (self.sfc * 1.1 * self.total_power 
		                      * mission.duration)
		self.empty_weight = sum(self.mass_estimation(
			power=power_msl).values())
		self.mtow = (self.empty_weight + self.fuel_mass 
		                 + mission.payload + mission.crew_mass)


	def mass_estimation(self, power):
		""" Estimate component masses of medium helicopters. [pp.408-416]

		Requires:
			Fuel mass (self.fuel_mass)
		"""
		# Wetted fuselage surface [m^2]
		wetted_surface = 59.09386 * np.exp(0.0000194463 * self.mtow)
		
		# Mass estimation [kg]
		m_mr = (33 * self.main_rotor.radius * self.main_rotor.chord 
		        * self.main_rotor.number_of_blades + 16)
		m_tr = 0.003942 * self.mtow + 5.66
		m_fus_t = 0.11907 * self.mtow - 66.666
		m_prop = 1.83 * (133.8 + 0.1156 * power * 1e-3)
		m_drive = (0.00000166 * (0.9 * self.mtow) ** 2 
		           + 0.087780096 * self.mtow - 113.81241656)
		m_tanks = 164.751 * np.log(self.fuel_mass / 2.948) - 751.33
		m_fcs = 95.6368 * np.exp(0.000111114 * self.mtow)
		m_i = 25.444 * np.log(power / 0.7457e3) - 141.62
		m_hyd = 0.003258 * self.mtow + 5.24
		m_el = 218.496 * np.log(wetted_surface / 0.092903) - 1267.49
		m_av = 113.4 + self.special_equipment
		m_furn = 0.854 * wetted_surface + 9.98 * self.number_of_seats - 4.54
		m_ac_ai = 55.542 * np.log(10.7369 * wetted_surface) - 331.21
		m_l_h = 38
		
		# Landing gear
		if self.landing_gear_type == 'Skids':
			m_lg = (0.011113004 * (0.9 * self.mtow / 0.453592) ** 0.8606 
			        * self.main_rotor.number_of_blades ** 0.8046)
		elif self.landing_gear_type == 'Rigid wheels':
			m_lg = (0.187333496 * (0.9 * self.mtow / 0.453592) ** 0.6662 
			        * self.number_of_legs ** 0.536)
		elif self.landing_gear_type == 'Retractable wheels':
			m_lg = (0.187333496 * (0.9 * self.mtow / 0.453592) ** 0.6662 
			        * 2 ** 0.1198 * self.number_of_legs ** 0.536)
		else: raise ValueError('Invalid landing gear type.')

		# Empty weight components [kg]
		return {
			'main rotor': m_mr, 'tail rotor': m_tr, 'fuselage and tail': m_fus_t, 
			'landing gear': m_lg, 'engines': m_prop, 'transmission': m_drive, 
			'fuel tanks': m_tanks, 'flight control systems': m_fcs, 
			'instruments': m_i, 'hydraulic systems': m_hyd, 
			'electrical systems': m_el, 'avionics': m_av, 'furnishing': m_furn,
			'air conditioning/anti-ice': m_ac_ai, 'loading and handling': m_l_h}


	def plot_powers(self, mission: Mission, max_velocity: float):
		""" Plot powers over a range of horizontal flight speeds for the current 
		MTOW. Level flight approximation is applied.
		"""
		# List of speeds (x-axis)
		speeds = np.linspace(0, max_velocity, 100)

		# Initialize parameters
		P_i, P_0, P_p, P_tr, P_tl, P_a, P = [], [], [], [], [], [], []
		weight = self.mtow * self.gravity
		mission.climb_angle = 0
		self.alpha = 0

		for V in speeds:

			# Flight state
			mission.flight_speed = V
			advance_ratio = V / self.main_rotor.tip_velocity
			induced_velocity = self.main_rotor.get_induced_velocity_level(
				mission.density, V, thrust=weight)

			# Power calculation [kW]
			P_i.append(self.main_rotor.kappa * weight * induced_velocity 
			           * 1e-3)
			P_0.append(self.main_rotor.get_profile_power(
				mission.density, advance_ratio) * 1e-3)
			P_p.append(self.get_parasite_power() * 1e-3)
			main_rotor_power = P_i[-1] + P_0[-1] + P_p[-1]
			P_tr.append(self.tail_rotor.power_fraction * main_rotor_power)
			P_tl.append(((1 / self.eta_transmission) - 1) * main_rotor_power)
			P_a.append(self.accessory_power * 1e-3)
			P.append(main_rotor_power + P_tr[-1] + P_tl[-1] + P_a[-1])

		# Plot
		fig, ax = plt.subplots(figsize=(8,5))
		ax.set_title('Power in level flight', fontweight='bold')
		ax.plot(speeds, P, label=r'$P$ Total power')
		ax.plot(speeds, P_i, label=r'$P_i$ Induced power')
		ax.plot(speeds, P_0, label=r'$P_0$ Profile power')
		ax.plot(speeds, P_p, label=r'$P_p$ Parasite power')
		ax.plot(speeds, P_tr, label=r'$P_{TR}$ Tail rotor power')
		ax.plot(speeds, P_tl, label=r'$P_{tl}$ Transmission losses')
		ax.plot(speeds, P_a, label=r'$P_a$ Accessory power')
		ax.set_xlabel('Velocity [m/s]')
		ax.set_ylabel('Power [kW]')
		fig.tight_layout()
		ax.legend()
		ax.grid()
		plt.show()

		# Reset mission segment
		mission.set_mission_segment(0)


	def plot_mtow_convergence(self, mtow_list):
		""" Plot MTOW over the iterations within the design loop.
		"""
		# Figure, plot
		fig, ax = plt.subplots(figsize=(8,5))
		ax.set_title('MTOW convergence', fontweight='bold')
		ax.plot(range(len(mtow_list)), mtow_list, 'k')
		ax.set_xlabel('Iterations')
		ax.set_ylabel('MTOW [kg]')
		fig.tight_layout()
		ax.grid()
		plt.show()


