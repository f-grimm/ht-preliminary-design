# -*- coding: utf-8 -*-
"""
Created on 2022-02-12

@author: Fabian Grimm (f.grimm@tum.de)
"""

import numpy as np

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
	helicopter configurations.

	Requires:
		Main rotor

	TODO: 
		Solidity: Blade loading diagram instead of constant chord.
		EW models [Prouty]
	"""
	# ---- Main rotor design [p.98-172] ---------------------------------------

	# Rotor radius for min. power in hover (more factors should be considered)
	rotor = aircraft.main_rotor
	rotor.solidity = rotor.get_solidity()
	weight = aircraft.mtow * aircraft.gravity
	rotor.radius = rotor.get_min_power_radius(
		aircraft.mission.density, thrust=weight)

	# ---- Preliminary performance estimation [p.175] -------------------------

	# Hover power as first estimate
	induced_power = rotor.get_induced_power_hover(
		aircraft.mission.density, thrust=weight)
	profile_power = rotor.get_profile_power(
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
	rotor.disc_loading = rotor.get_disc_loading(thrust=weight)
	rotor.figure_of_merit = ((induced_power / rotor.kappa) 
	                         / aircraft.hover_power)


def iterate_second_sizing(aircraft: object):
	""" Perform one step in the second sizing loop for conventional 
	helicopter configurations.

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
	tr = aircraft.tail_rotor
	tr.radius = 0.4 * np.sqrt(2.2 * aircraft.mtow / 1000)
	tr.solidity = tr.get_solidity()

	# ---- Refined performance calculation [pp. 251 - 321] --------------------

	# Flight state calculation
	drag = aircraft.get_fuselage_drag()
	aircraft.alpha = aircraft.get_angle_of_attack(drag)
	rotor = aircraft.main_rotor
	aircraft.advance_ratio = aircraft.mission.flight_speed / rotor.tip_velocity
	thrust = aircraft.get_thrust(drag)
	induced_velocity = rotor.get_induced_velocity(
		aircraft.mission.density, aircraft.mission.flight_speed, 
		aircraft.alpha, thrust)

	# Power calculation
	induced_power = rotor.kappa * thrust * induced_velocity
	profile_power = rotor.get_profile_power(
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
	aircraft.mtow = get_initial_mtow(aircraft)

	if logs:
		print(f'Mission segment: 1')
		print(f'Height: {aircraft.mission.height:.2f} m\n-----')
		print(f'First MTOW estimation: {aircraft.mtow:10.2f} kg')

	iterate_first_sizing(aircraft)
	iterate_second_sizing(aircraft)

	if logs:
		print(f'Hover power: {(aircraft.hover_power / 1000):20.2f} kW')
		print(f'Tail rotor radius: {aircraft.tail_rotor.radius:14.2f} m')
