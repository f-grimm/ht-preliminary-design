# -*- coding: utf-8 -*-
"""
Created on 2022-02-12

@author: Fabian Grimm (f.grimm@tum.de)
"""

import numpy as np
import matplotlib.pyplot as plt

"""
"""

def preliminary_design(aircraft: object, mission_segment: int, logs=True):
	""" Preliminary design of conventional helicopter configurations.

	TODO:
		Design loop
	"""
	aircraft.mission.set_mission_segment(mission_segment)

	if logs:
		print('\nPreliminary Design\n' + '-' * 18)
		print(f'\nAicraft: {aircraft.name}')
		print(f'Mission: {aircraft.mission.name}')
		print(f'Segment: {aircraft.mission.segment} '
		      + f'({aircraft.mission.height:.0f}m, '
		      + f'{aircraft.mission.duration:.2f}h, '
		      + f'{aircraft.mission.flight_speed:.0f}m/s)')
		print(f'\nEstimate initial MTOW'
		      + f'\n(SFC {aircraft.sfc * 1e3:.2f}, '
		      + f'EWR {aircraft.empty_weight_ratio:.2f})')

	aircraft.mtow = get_initial_mtow(aircraft)
		
	if logs:
		print(f'\n - MTOW: {aircraft.mtow:27.2f} kg')
		print('\nIterate first sizing'
		      + '\n(Hover, min. power)')

	iterate_first_sizing(aircraft)

	if logs:
		print(f'\n - Main rotor radius: {aircraft.main_rotor.radius:14.2f} m')
		print(f' - Hover power: {(aircraft.hover_power * 1e-3):20.2f} kW')
		print(f' - MTOW: {aircraft.mtow:27.2f} kg')
		print('\nIterate second sizing'
		      + '\n(Refined performance)')

	iterate_second_sizing(aircraft)

	if logs:
		print(f'\n - Tail rotor radius: {aircraft.tail_rotor.radius:14.2f} m')
		print(f' - Drag: {aircraft.drag:27.2f} N')
		print(f' - Thrust: {aircraft.thrust:25.2f} N')
		print(f' - Angle of attack: {aircraft.alpha * 180 / np.pi:16.2f} deg')
		print(f' - Advance ratio: {aircraft.advance_ratio:18.2f} -')
		print(' - Induced velocity: '
		      + f'{aircraft.main_rotor.induced_velocity:15.2f} m/s')
		print(f' - MTOW: {aircraft.mtow:27.2f} kg')
		print('')


def get_initial_mtow(aircraft: object):
	""" Estimate the initial maximum take-off weight based on the mission 
	profile; (SFC and empty weight ratio are constant). [p.90-91]
	"""
	fuel_mass = (aircraft.sfc * aircraft.power_available 
	             * aircraft.mission.duration)
	useful_load = (fuel_mass + aircraft.mission.payload 
	               + aircraft.mission.crew_mass)

	# MTOW [kg]
	return useful_load / (1 - aircraft.empty_weight_ratio)


def iterate_first_sizing(aircraft: object):
	""" Perform one step in the first sizing loop [p. 86]

	TODO: 
		Solidity: Blade loading diagram instead of constant chord.
	"""
	# ---- Main rotor sizing [p.98-172] ---------------------------------------

	# Rotor radius for min. power in hover (more factors should be considered)
	aircraft.main_rotor.solidity = aircraft.main_rotor.get_solidity()
	weight = aircraft.mtow * aircraft.gravity
	aircraft.main_rotor.radius = aircraft.main_rotor.get_min_power_radius(
		aircraft.mission.density, thrust=weight)

	# ---- Preliminary performance estimation [p.175] -------------------------

	# Hover power
	induced_power = aircraft.main_rotor.get_induced_power_hover(
		aircraft.mission.density, thrust=weight)
	profile_power = aircraft.main_rotor.get_profile_power(
		aircraft.mission.density, advance_ratio=0)
	aircraft.hover_power = induced_power + profile_power

	# ---- Mass estimation incl. fuel [p.177-179] -----------------------------

	# Fuel consumption in hover, empty weight
	aircraft.fuel_mass = (aircraft.sfc * 1.1 * aircraft.hover_power 
	                      * aircraft.mission.duration)
	aircraft.empty_weight = sum(mass_estimation(
		aircraft, power=aircraft.hover_power).values())
	aircraft.mtow = (aircraft.empty_weight + aircraft.fuel_mass 
	                 + aircraft.mission.payload + aircraft.mission.crew_mass)

	# Disc loading and FM [p.181-183]
	weight = aircraft.mtow * aircraft.gravity
	aircraft.main_rotor.disc_loading = aircraft.main_rotor.get_disc_loading(
		thrust=weight)
	aircraft.main_rotor.figure_of_merit = (
		(induced_power / aircraft.main_rotor.kappa) / aircraft.hover_power)


def iterate_second_sizing(aircraft: object):
	""" Perform one step in the second sizing loop. [p. 86]

	TODO: 
		Solidity: Empirical reference instead of constant chord.
		How should the download factor be considered? What counts as slow
			flight?
		Engine calc., SFC = f(P) [p. 310]
		Conditions at the beginning of each segment assumed. Half altitude for 
			climb?
	"""
	# ---- Tail rotor design [p.204-213] --------------------------------------

	# Tail rotor design according to Layton
	aircraft.tail_rotor.radius = 0.4 * np.sqrt(2.2 * aircraft.mtow * 1e-3)
	aircraft.tail_rotor.solidity = aircraft.tail_rotor.get_solidity()

	# ---- Refined performance calculation [pp. 251 - 321] --------------------

	# Flight state
	aircraft.drag = aircraft.get_fuselage_drag()
	aircraft.alpha = aircraft.get_angle_of_attack()
	aircraft.thrust = aircraft.get_thrust()
	aircraft.advance_ratio = (aircraft.mission.flight_speed 
	                          / aircraft.main_rotor.tip_velocity)
	aircraft.main_rotor.induced_velocity = (
		aircraft.main_rotor.get_induced_velocity(
			aircraft.mission.density, aircraft.mission.flight_speed, 
			aircraft.alpha, aircraft.thrust))

	# Power calculation
	induced_power = (aircraft.main_rotor.kappa * aircraft.thrust 
	                 * aircraft.main_rotor.induced_velocity)
	profile_power = aircraft.main_rotor.get_profile_power(
		aircraft.mission.density, aircraft.advance_ratio)
	parasite_power = aircraft.get_parasite_power()
	climb_power = aircraft.get_climb_power()
	main_rotor_power = (induced_power + profile_power + parasite_power 
	                    + climb_power)
	tail_rotor_power = aircraft.tail_rotor.power_fraction * main_rotor_power
	transmission_losses = (((1 / aircraft.eta_transmission) - 1) 
	                       * main_rotor_power)
	aircraft.total_power = (main_rotor_power + tail_rotor_power 
	                        + transmission_losses + aircraft.accessory_power)

	# Engine requirement
	temperature_ratio = aircraft.mission.temperature / 288.15
	pressure_ratio = aircraft.mission.pressure / 101325
	power_msl = (aircraft.total_power * np.sqrt(temperature_ratio) 
	             / pressure_ratio)

	# ---- Refined mass estimation [p.324-325] --------------------------------
	
	# Fuel mass, empty weight
	aircraft.fuel_mass = (aircraft.sfc * 1.1 * aircraft.total_power 
	                      * aircraft.mission.duration)
	aircraft.empty_weight = sum(mass_estimation(
		aircraft, power=power_msl).values())
	aircraft.mtow = (aircraft.empty_weight + aircraft.fuel_mass 
	                 + aircraft.mission.payload + aircraft.mission.crew_mass)


def mass_estimation(aircraft: object, power):
	""" Estimate component masses of medium helicopters. [pp. 408-416]

	Requires:
		Fuel mass (aircraft.fuel_mass)
	"""
	# Wetted fuselage surface [m^2]
	wetted_surface = 59.09386 * np.exp(0.0000194463 * aircraft.mtow)
	
	# Mass estimation [kg]
	m_mr = (33 * aircraft.main_rotor.radius * aircraft.main_rotor.chord 
	        * aircraft.main_rotor.number_of_blades + 16)
	m_tr = 0.003942 * aircraft.mtow + 5.66
	m_fus_t = 0.11907 * aircraft.mtow - 66.666
	m_prop = 1.83 * (133.8 + 0.1156 * power * 1e-3)
	m_drive = (0.00000166 * (0.9 * aircraft.mtow) ** 2 
	           + 0.087780096 * aircraft.mtow - 113.81241656)
	m_tanks = 164.751 * np.log(aircraft.fuel_mass / 2.948) - 751.33
	m_fcs = 95.6368 * np.exp(0.000111114 * aircraft.mtow)
	m_i = 25.444 * np.log(power / 0.7457e3) - 141.62
	m_hyd = 0.003258 * aircraft.mtow + 5.24
	m_el = 218.496 * np.log(wetted_surface / 0.092903) - 1267.49
	m_av = 113.4 + aircraft.special_equipment
	m_furn = 0.854 * wetted_surface + 9.98 * aircraft.number_of_seats - 4.54
	m_ac_ai = 55.542 * np.log(10.7369 * wetted_surface) - 331.21
	m_l_h = 38
	
	# Landing gear (default: skids)
	if aircraft.landing_gear_type == 'Skids':
		m_lg = (0.011113004 * (0.9 * aircraft.mtow / 0.453592) ** 0.8606 
		        * aircraft.main_rotor.number_of_blades ** 0.8046)
	elif aircraft.landing_gear_type == 'Rigid wheels':
		m_lg = (0.187333496 * (0.9 * aircraft.mtow / 0.453592) ** 0.6662 
		        * aircraft.number_of_legs ** 0.536)
	elif aircraft.landing_gear_type == 'Retractable wheels':
		m_lg = (0.187333496 * (0.9 * aircraft.mtow / 0.453592) ** 0.6662 
		        * 2 ** 0.1198 * n_legs ** 0.536)
	else: raise ValueError('Invalid landing gear type.')

	# Empty weight components [kg]
	return {
		'main rotor': m_mr, 'tail rotor': m_tr, 'fuselage and tail': m_fus_t, 
		'landing gear': m_lg, 'engines': m_prop, 'transmission': m_drive, 
		'fuel tanks': m_tanks, 'flight control systems': m_fcs, 
		'instruments': m_i, 'hydraulic systems': m_hyd, 
		'electrical systems': m_el, 'avionics': m_av, 'furnishing': m_furn,
		'air conditioning/anti-ice': m_ac_ai, 'loading and handling': m_l_h}


def plot_powers(aircraft: object, max_velocity):
	""" Plot powers over a range of horizontal flight speeds for the current 
	MTOW. Level flight approximation is applied.
	"""
	# List of speeds (x-axis)
	speeds = np.linspace(0, max_velocity, 100)

	# Initialize parameters
	P_i, P_0, P_p, P_tr, P_tl, P_a, P = [], [], [], [], [], [], []
	weight = aircraft.mtow * aircraft.gravity
	aircraft.mission.climb_angle = 0
	aircraft.alpha = 0

	for V in speeds:

		# Flight state
		aircraft.mission.flight_speed = V
		advance_ratio = V / aircraft.main_rotor.tip_velocity
		induced_velocity = aircraft.main_rotor.get_induced_velocity_level(
			aircraft.mission.density, V, thrust=weight)

		# Power calculation [kW]
		P_i.append(aircraft.main_rotor.kappa * weight * induced_velocity 
		           * 1e-3)
		P_0.append(aircraft.main_rotor.get_profile_power(
			aircraft.mission.density, advance_ratio) * 1e-3)
		P_p.append(aircraft.get_parasite_power() * 1e-3)
		main_rotor_power = P_i[-1] + P_0[-1] + P_p[-1]
		P_tr.append(aircraft.tail_rotor.power_fraction * main_rotor_power)
		P_tl.append(((1 / aircraft.eta_transmission) - 1) * main_rotor_power)
		P_a.append(aircraft.accessory_power * 1e-3)
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
	aircraft.mission.set_mission_segment(0)
