# -*- coding: utf-8 -*-
"""
Created on 2022-02-12

@author: Fabian Grimm (f.grimm@tum.de)
"""

import numpy as np
import matplotlib.pyplot as plt

"""
"""

def get_initial_mtow(aircraft: object):
	""" Estimate the initial maximum take-off weight based on the mission 
	profile. [p.90-91]
	"""
	fuel_mass = (aircraft.sfc * aircraft.power_available 
	             * aircraft.mission.duration)
	useful_load = (fuel_mass + aircraft.mission.payload 
	               + aircraft.mission.crew_mass)

	# MTOW [kg]
	return useful_load / (1 - aircraft.empty_weight_ratio)


def iterate_first_sizing(aircraft: object):
	""" Perform one step in the first sizing loop for conventional 
	helicopter configurations. [p. 86]

	Requires:
		Main rotor

	TODO: 
		Solidity: Blade loading diagram instead of constant chord.
		EW models [Prouty]
	"""
	# ---- Main rotor design [p.98-172] ---------------------------------------

	# Rotor radius for min. power in hover (more factors should be considered)
	aircraft.main_rotor.solidity = aircraft.main_rotor.get_solidity()
	weight = aircraft.mtow * aircraft.gravity
	aircraft.main_rotor.radius = aircraft.main_rotor.get_min_power_radius(
		aircraft.mission.density, thrust=weight)

	# ---- Preliminary performance estimation [p.175] -------------------------

	# Hover power as first estimate
	induced_power = aircraft.main_rotor.get_induced_power_hover(
		aircraft.mission.density, thrust=weight)
	profile_power = aircraft.main_rotor.get_profile_power(
		aircraft.mission.density, advance_ratio=0)
	aircraft.hover_power = induced_power + profile_power

	# ---- Mass estimation incl. fuel [p.177-179] -----------------------------

	# Empty weight estimation, fuel demand in hover
	empty_weight = aircraft.empty_weight_ratio * aircraft.mtow
	fuel_mass = (aircraft.sfc * 1.1 * aircraft.hover_power 
	             * aircraft.mission.duration)
	aircraft.mtow = (empty_weight + fuel_mass + aircraft.mission.payload 
	                 + aircraft.mission.crew_mass)

	# Disc loading and FM [p.181-183]
	weight = aircraft.mtow * aircraft.gravity
	aircraft.main_rotor.disc_loading = aircraft.main_rotor.get_disc_loading(
		thrust=weight)
	aircraft.main_rotor.figure_of_merit = (
		(induced_power / aircraft.main_rotor.kappa) / aircraft.hover_power)


def iterate_second_sizing(aircraft: object):
	""" Perform one step in the second sizing loop for conventional 
	helicopter configurations. [p. 86]

	Requires:
		Main rotor
		Tail rotor

	TODO: 
		Solidity: Empirical reference instead of constant chord.
		How should the download factor be considered? What counts as slow
			flight?
	"""
	# ---- Tail rotor design [p.204-213] --------------------------------------

	# Tail rotor design according to Layton
	aircraft.tail_rotor.radius = 0.4 * np.sqrt(2.2 * aircraft.mtow / 1000)
	aircraft.tail_rotor.solidity = aircraft.tail_rotor.get_solidity()

	# ---- Refined performance calculation [pp. 251 - 321] --------------------

	# Flight state
	aircraft.drag = aircraft.get_fuselage_drag()
	aircraft.alpha = aircraft.get_angle_of_attack()
	aircraft.thrust = aircraft.get_thrust()
	aircraft.advance_ratio = (aircraft.mission.flight_speed 
	                          / aircraft.main_rotor.tip_velocity)
	induced_velocity = aircraft.main_rotor.get_induced_velocity(
		aircraft.mission.density, aircraft.mission.flight_speed, 
		aircraft.alpha, aircraft.thrust)

	# Power calculation
	induced_power = (aircraft.main_rotor.kappa * aircraft.thrust 
	                 * induced_velocity)
	profile_power = aircraft.main_rotor.get_profile_power(
		aircraft.mission.density, aircraft.advance_ratio)
	parasite_power = aircraft.get_parasite_power()
	climb_power = aircraft.get_climb_power()
	main_rotor_power = (induced_power + profile_power + parasite_power 
	                    + climb_power)
	tail_rotor_power = 0.05 * main_rotor_power
	transmission_losses = (((1 / aircraft.eta_transmission) - 1) 
	                       * main_rotor_power)
	total_power = (main_rotor_power + tail_rotor_power + transmission_losses 
	               + aircraft.accessory_power)

	# ...
	
	# Selection of engine and gearbox
	# Mass estimation incl. fuel


def preliminary_design(aircraft: object, logs=True):
	""" Preliminary design of conventional helicopter configurations.

	Requires:
		Main rotor
		Tail rotor

	TODO:
		Design loop
	"""
	if logs:
		print('\nPreliminary Design\n' + '-' * 18)
		print(f'\nMission segment: {aircraft.mission.segment} '
		      + f'({aircraft.mission.height:.2f}m)')
		print(f'\nEstimate initial MTOW'
		      + f'\n(SFC {aircraft.sfc * 1e3:.2f}, '
		      + f'{aircraft.mission.duration:.2f}h, '
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
		print('')


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

		# Power calculation
		P_i.append(aircraft.main_rotor.kappa * weight * induced_velocity 
		           * 1e-3)
		P_0.append(aircraft.main_rotor.get_profile_power(
			aircraft.mission.density, advance_ratio) * 1e-3)
		P_p.append(aircraft.get_parasite_power() * 1e-3)
		main_rotor_power = P_i[-1] + P_0[-1] + P_p[-1]
		P_tr.append(0.05 * main_rotor_power)
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
