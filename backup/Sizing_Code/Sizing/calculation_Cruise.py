#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 15:18:19 2020

@author: christian
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

def power_induced(fstate, comp):
    ct_sig = comp['ct_sig']
    sigma = comp['sigma']
    radius = comp['radius']
    omega = fstate['omega']
    thrust = thrust_from_omega(fstate, comp)
    lamb_hov = math.sqrt(ct_sig * sigma / 2)
    lamb_ind = optimize.newton(inflow_equation, lamb_hov, args=(fstate, comp))    
    return lamb_ind * omega * radius * thrust

def inflow_equation(lamb_ind, fstate, comp):
    radius = comp['radius']
    ct_sig = comp['ct_sig']
    sigma = comp['sigma']
    omega = fstate['omega']
    V = fstate['V']
    mu = V / (omega*radius)
    mu_z = 0.
    lamb_hov = math.sqrt(ct_sig * sigma/ 2)
    return lamb_ind * math.sqrt( mu**2 + (mu_z + lamb_ind)**2 ) - lamb_hov**2

def power_profile(fstate, comp):
    radius = comp['radius']
    c_d0 = comp['c_d0']
    sigma = comp['sigma']
    omega = fstate['omega']
    V = fstate['V']
    dens = fstate['dens']
    
    A = radius**2*math.pi
    mu = V / (omega*radius)
    p0h = sigma/8*c_d0*(dens*A*(omega*radius)**3)
    p0 = p0h * (1 + 4.65*mu**2)
    return p0

def thrust_from_omega(fstate, comp):
    omega = fstate['omega'] 
    dens = fstate['dens']
    radius = comp['radius']
    sigma = comp['sigma']
    ct_sig = comp['ct_sig']
    A = math.pi * radius**2
    thrust = ct_sig * sigma * dens * A * (omega*radius)**2
    return thrust
def calcSigma(comp):
    n_bld = comp['n_bld']
    chord = comp['chord']
    radius = comp['radius']
    return n_bld*chord*radius/(math.pi * radius**2)
def calcLift(fstate, comp):
    dens = fstate['dens']
    V = fstate['V']
    c_l = fstate['c_l']
    chord = comp['chord']
    l = comp['l']
    S_ref = chord*l

    return  0.5 * dens * S_ref * V**2 * c_l

def calcDrag(fstate, comp):
    dens = fstate['dens']
    V = fstate['V']
    c_l = fstate['c_l']
    if comp['type'] == 'wing':
        chord = comp['chord']
        l = comp['l']
        S_ref = chord*l
        e = comp['e']
        k = chord/(l*math.pi*e)
        return 0.5 * dens * S_ref * V**2 * k * c_l**2
    elif comp['type'] == 'fuselage':
        efpa = comp['efpa']
        return 0.5 * dens * V**2 * efpa
    else:
        print("Error")
        return 0.
def calcGravity(m):
    return m * 9.81

def get_matrix_g_to_cg(theta, phi):
    """
    returns the transformation matrix from the erdlotfestes system to the
    body system.
    Parameters
    ----------
    theta : pitch angle
    phi : bank angle
     
    Returns
    -------
    trafomatrix g to cg
    """
    return np.array([[ math.cos(theta),     math.sin(theta)*math.cos(phi), -1.*math.sin(theta)*math.cos(phi) ], 
                     [              0.,                   math.cos(theta),                     math.sin(phi) ], 
                     [ math.sin(theta), -1.*math.cos(theta)*math.sin(phi),     math.cos(theta)*math.cos(phi) ]])

def sumForcesCruise(sol, m, fstate, components, results):

    powers_ind = []
    powers_prf = []
    powers_par = []
    
    drags_ind = []
    drags_par = []

    lifts = []
    drags = []
    thrusts = []
    gravities = []

    atmo = calcAtmo(fstate['h'])
    fstate['temp'] = atmo[0]
    fstate['dens'] = atmo[1]
    fstate['a'] = atmo[2]
    
    fstate['omega'] = sol[0]
    fstate['c_l'] = sol[1]
    
    for comp in components:
        if comp['type'] == 'propeller':
            comp['sigma'] = calcSigma(comp)
            # --- calculate generated thrust
            thrusts.append(thrust_from_omega(fstate, comp))
            # --- sum force and moment vectors and append to results list
            powers_ind.append(power_induced(fstate, comp))
            powers_prf.append(power_profile(fstate, comp))
            results['v_tip'] = comp['radius'] * fstate['omega']
        elif comp['type'] == 'wing':
            lifts.append(calcLift(fstate, comp))
            drags.append(calcDrag(fstate, comp))
            powers_ind.append(drags[-1] * fstate['V'])
            drags_ind.append(drags[-1])
            
        elif comp['type'] == 'fuselage':
            gravities.append(calcGravity(m))
            drags.append(calcDrag(fstate, comp))
            powers_par.append(drags[-1] * fstate['V'])
            drags_par.append(drags[-1])
            
            
    results['lift'] = sum(lifts)
    results['drag'] = sum(drags)
    results['thrust'] = sum(thrusts)
    results['gravity'] = sum(gravities)
    
    thrust_b = np.array([1., 0., 0.]) * results['thrust']
    lift_b = np.array([0., 0., -1.]) * results['lift']
    drag_b = np.array([-1., 0., 0.]) * results['drag']
    gravity_b = get_matrix_g_to_cg(fstate['gamma'], 0) @ np.array([0., 0., 1.]) * results['gravity']
            
    results['power_ind'] = sum(powers_ind)
    results['power_prf'] = sum(powers_prf)
    results['power_par'] = sum(powers_par)
    results['power_ges'] = results['power_ind'] + results['power_prf'] + results['power_par']
    results['drag_ind'] = sum(drags_ind)
    results['drag_par'] = sum(drags_par)
    
    results['h_dot'] = fstate['V'] * math.sin(fstate['gamma'])
    results['atmo'] = atmo
    
    residual = thrust_b + lift_b + drag_b + gravity_b

    return residual

def getPowerCruise(m, dData):
    # set x0 values
    omega = 200.           # --- angular velocity of all props
    c_l = 0.5            # --- Lift coefficient
    mu = 0.                # --- bank angle
    
    # Initialize vector (list) for x0
    x0 = [omega, c_l, mu]
    components = dData['configuration']['cruise']
    V = dData['mission']['cruiseVelo']
    h = dData['mission']['altitude']
    results = {}
    flightstate = {'V': V, 'gamma': 0., 'h': h}
    sol = optimize.root(sumForcesCruise, x0, args=(m, flightstate, components, results))
    dData['results']['powers']['cruisePower'] = results['power_ges']
    dData['results']['ltod'] = results['lift']/results['drag']
    dData['results']['c_l'] = sol.x[1]
    return dData

def stationaryCrs(dData):
    components = dData['configuration']['cruise']
    m = dData['masses']['dgm']
    
    # set x0 values
    omega = 200.           # --- angular velocity of all props
    c_l = 0.5            # --- Lift coefficient
    mu = 0.                # --- bank angle
    
    x0 = [omega, c_l, mu]
    
    # Define flightstate
    
    V = 65.
    gamma = 0./180*math.pi
    h = 1000.
    
    var = 'V'
    
    if var == 'gamma':
        values = np.arange(0,50/180*math.pi,0.01)
    elif var == 'V':
        values = np.arange(30,71, 1.)
    elif var == 'h':
        values = np.arange(0,10000, 100)
    
    omegas = []
    c_ls = []
    pInd = []
    pPrf = []
    pPar = []
    pGes = []
    vTips = []
    drags_ges = []
    drags_ind = []
    drags_par = []
    lifts = []
    ltodrags = []
    
    for value in values:
        flightstate = {'type': 'cruise', 'V': V, 'gamma': gamma, 'h': h}
        flightstate.update({var: value})
        results = {}
        sol = optimize.root(sumForcesCruise, x0, args=(m, flightstate, components, results))
        a = results['atmo'][2]
        omegas.append(sol.x[0])
        c_ls.append(sol.x[1])
        
        pInd.append(results['power_ind'] / 1000.)
        pPrf.append(results['power_prf'] / 1000.)
        pPar.append(results['power_par'] / 1000.)
        pGes.append(pInd[-1] + pPrf[-1] + pPar[-1])
        
        vTips.append(results['v_tip']/a) 
        
        lifts.append(results['lift'] / 1000.)
        drags_ges.append(results['drag'] / 1000.)
        drags_ind.append(results['drag_ind'] / 1000.)
        drags_par.append(results['drag_par'] / 1000.)
        ltodrags.append(lifts[-1]/drags_ges[-1])
        
    powers = {'power_ges': pGes, 'power_ind': pInd, 'power_prf': pPrf, 'power_par': pPar} 
    forces = {'drags_ges': drags_ges, 'drags_ind': drags_ind, 'drags_par': drags_par,
              'ltodrags': ltodrags}
        
    results = {'omegas': omegas, 'c_ls': c_ls, 'powers': powers, 'forces': forces}
       
    generatePlot(flightstate, var, values, dData, results)
    
    return dData