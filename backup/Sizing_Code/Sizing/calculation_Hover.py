#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 11:05:13 2019

@author: ge56zeg
"""
import math
import numpy as np
from scipy import optimize
from plotting import generatePlot

def calcAtmo(h):
    """
    Calculates the atmospheric temperature, the density and the speed of sound
    at a given altitude h in the Troposphere.

    Parameters
    ----------
    h : altitude in m above NN

    Returns
    -------
    T : temperature in K
    rho : density in kg/m^3
    a : speed of sound in m/s

    """
    gamma = -0.0065     
    TMSL = 288.15       #Temperatur auf MSL in K
    n = 1.235           
    rhoMSL = 1.225      #Luftdichte auf MSL in kg/m^3
    rE = 6356766        #Erdradius in m
    R = 287.05          #Gaskonstante für Luft
    kappa = 1.402       #Isentropenkoeffizient für Luft
    
    HG = rE * h / (rE + h)
    temp = TMSL + gamma * HG
    rho = rhoMSL * (temp / TMSL)**(1/(n-1))
    a = math.sqrt(kappa*R*temp)
    return temp, rho, a

def thrust_from_omega(fstate, comp):
    omega = fstate['omega']
    radius = comp['radius']
    dens = fstate['dens']
    A = math.pi * radius**2
    sigma = comp['sigma']
    thrust = comp['ct_sig'] * sigma * dens * A * (omega*radius)**2
    return thrust

def calc_weight(mass):
    return mass * 9.80665

def power_induced(fstate, comp):
    return vel_induced(fstate, comp) * thrust_from_omega(fstate, comp)

def vel_induced(fstate, comp):
    vh = vel_induced_hover(fstate, comp)
    V = fstate['V']
    return -V/2+math.sqrt((V/2)**2+vh**2)

def vel_induced_hover(fstate, comp):
    radius = comp['radius']
    dens = fstate['dens']
    return math.sqrt(thrust_from_omega(fstate, comp) / (2*dens*math.pi*radius**2))

def power_climb(fstate, comp):
    return thrust_from_omega(fstate, comp) * fstate['V']

def power_profile(fstate, comp):
    omega = fstate['omega']
    dens = fstate['dens']
    radius = comp['radius']
    A = math.pi * radius**2
    sigma = comp['sigma']
    p0h = sigma/8*comp['c_d0']*(dens*A*(omega*radius)**3)
    p0 = p0h
    return p0

def get_matrix_a_to_cg(alpha):
    return np.array([[math.cos(alpha), 0., -1.*math.sin(alpha)],
                     [0.,              1.,                  0.],
                     [math.sin(alpha), 0.,     math.cos(alpha)]])

def get_matrix_g_to_cg(theta, phi):
    return np.array([[ math.cos(theta),     math.sin(theta)*math.cos(phi), -1.*math.sin(theta)*math.cos(phi) ], 
                     [              0.,                   math.cos(theta),                     math.sin(phi) ], 
                     [ math.sin(theta), -1.*math.cos(theta)*math.sin(phi),     math.cos(theta)*math.cos(phi) ]])
   
def get_angles_aer_cg(vel_cg):
    vel_abs = math.sqrt(vel_cg[0]**2 + vel_cg[1]**2 + vel_cg[2]**2)
    beta = math.asin(vel_cg[1] / vel_abs)
    alpha = -1*math.acos(vel_cg[0] / math.cos(beta) / vel_abs)
    return alpha, beta
    
def inflow_equation(lamb_ind, omega, fstate, comp):
    radius = comp['radius']
    vel_cg = fstate['vel_cg']
    vel_abs = fstate['vel_abs']
    alpha = get_angles_aer_cg(vel_cg)[0]
    mu = vel_abs*math.cos(alpha) / (omega*radius)
    mu_z = -vel_abs*math.sin(alpha) / (omega*radius)
    lamb_hov = math.sqrt(comp['ct_sig'] / 2)
    return lamb_ind * math.sqrt( mu**2 + (mu_z + lamb_ind)**2 ) - lamb_hov**2

def calcSigma(comp):
    n_bld = comp['n_bld']
    chord = comp['chord']
    radius = comp['radius']
    return n_bld*chord*radius/(math.pi * radius**2)
    

def sumForcesHover(sol, m, fstate, components, results):
    
    forces = []
    powers_ind = []
    powers_clb = []
    powers_prf = []
    
    fstate['omega'] = sol[0]
    atmo = calcAtmo(fstate['h'])
    fstate['temp'] = atmo[0]
    fstate['dens'] = atmo[1]
    fstate['a'] = atmo[2]
    
    for comp in components:
        if comp['type'] == 'propeller':
            comp['sigma'] = calcSigma(comp)
            powerInd = power_induced(fstate, comp)
            powerClimb = power_climb(fstate, comp)
            powerPrf = power_profile(fstate, comp)

            powers_ind.append(powerInd)
            powers_prf.append(powerPrf)
            powers_clb.append(powerClimb)
            thrust = thrust_from_omega(fstate, comp)
            forces.append(-thrust)
        elif comp['type'] == 'fuselage':
            # --- set gravity vector and transform to cg system
            G = calc_weight(m)
            
            # --- sum force and moment vectors and append to results list 
            forces.append(G)
    
    results['power_ind'] = sum(powers_ind)
    results['power_prf'] = sum(powers_prf)
    results['power_clb'] = sum(powers_clb)
    results['power_ges'] = results['power_ind'] + results['power_prf'] + results['power_clb']
    residual = sum(forces)
    
    return residual

def getPowerHover(m, dData):
    # set x0 values
    omega = 200.      # rad/s

    # Initialize vector (list) for x0
    
    x0 = omega
    results = {}
    h = dData['mission']['altitude'] / 2
    V = dData['mission']['climbVelo']
    components = dData['configuration']['hover']
    
    flightstate = {'V': V, 'h': h}
    optimize.root(sumForcesHover, x0, args=(m, flightstate, components, results))
    dData['results']['powers']['climbPower'] = results['power_ges']
    
    flightstate = {'V': -V*1.5, 'h': h}
    optimize.root(sumForcesHover, x0, args=(m, flightstate, components, results))
    dData['results']['powers']['sinkPower'] = results['power_ges']
    return dData

def stationaryHov(dData):
    
    omega = 200.      # rad/s
    
    x0 = omega
    
    components = dData['configuration']['hover']
    # Define flightstate
    flightstate = {'type': 'hover', 'V': 5, 'h': 250}
    m = dData['masses']['dgm']
    
    velocities = range(-15, 16)
    
    pInd = []
    pClb = []
    pPrf = []
    pGes = []
    for vel in velocities:
        flightstate.update({'V': vel})
        results = {}
        optimize.root(sumForcesHover, x0, args=(m, flightstate, components, results))
        pInd.append(results['power_ind'] / 1000.)
        pPrf.append(results['power_prf'] / 1000.)
        pClb.append(results['power_clb'] / 1000.)
        pGes.append(results['power_ges'] / 1000.)
        
    powers = {'power_ind': pInd, 'power_clb': pClb, 'power_prf': pPrf, 'power_ges': pGes}
        
    results = {'powers': powers}
    
    
       
    generatePlot(flightstate, 'velocity', velocities, dData, results)
    return dData

