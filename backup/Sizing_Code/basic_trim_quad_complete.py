#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 11:05:13 2019

@author: gu95pug
"""
import sys
import math
import numpy as np

import matplotlib.pyplot as plt
from scipy import optimize

def thrust_from_omega(omega, dictFstate, dictRotor):##
    thrust=dictFstate['dens']*dictRotor['ct_sig']*dictRotor['n_bld']*dictRotor['chord']*omega**2*dictRotor['radius']**3##
    return thrust ### ### Tip: Suche nach einem Zusammenhang von Schub und Drehzahl

def calc_weight(mass):
    return mass * 9.80665

def calc_drag_aer(dictFstate,body):##
    v=math.sqrt(dictFstate['vel_x_g']**2+dictFstate['vel_y_g']**2+dictFstate['vel_z_g']**2)##
    drag_aer=0.5*dictFstate['dens']*v**2*body['efpa']##
    return drag_aer ### hhfhfh

def power_induced_hover(omega, dictFstate, dictRotor):
     return vel_induced(omega, dictFstate, dictRotor) * thrust_from_omega(omega, dictFstate, dictRotor)


def power_induced(omega, dictFstate, dictRotor):##
    ct=dictRotor['ct_sig']*dictRotor['n_bld']*dictRotor['chord']/(math.pi*dictRotor['radius'])##
    lamb_hov =math.sqrt(ct/2)##
    lamb_ind = optimize.newton(inflow_equation, lamb_hov, args=(omega, dictFstate, dictRotor))   
    return lamb_ind * omega * dictRotor['radius'] * thrust_from_omega(omega, dictFstate, dictRotor)

def power_profile(omega, dictFstate,dictRotor):##
    sigma=dictRotor['n_bld']*dictRotor['chord']/(math.pi*dictRotor['radius'])##
    A=math.pi*dictRotor['radius']**2##
    p_profile=sigma/8*dictRotor['c_d0']*(dictFstate['dens']*A*(omega*dictRotor['radius'])**3)##
    mu=calc_mu(omega, dictFstate, dictRotor)
    return p_profile * (1. + 4.65*mu**2)
                        ###
def vel_induced(omega, dictFstate, dictRotor):
    radius = dictRotor['radius']
    dens = dictFstate['dens']
    return math.sqrt(thrust_from_omega(omega, dictFstate, dictRotor) / (2*dens*math.pi*radius**2) )

# =============================================================================
# def get_matrix_a_to_cg(alpha):
#     return np.array([[math.cos(alpha), 0., -1.*math.sin(alpha)],
#                      [0.,              1.,                  0.],
#                      [math.sin(alpha), 0.,     math.cos(alpha)]])
# =============================================================================

def get_matrix_g_to_cg(theta, phi):
    return np.array([[ math.cos(theta),     math.sin(theta)*math.cos(phi), -1.*math.sin(theta)*math.cos(phi) ], 
                     [              0.,                   math.cos(theta),                     math.sin(phi) ], 
                     [ math.sin(theta), -1.*math.cos(theta)*math.sin(phi),     math.cos(theta)*math.cos(phi) ]])

def get_matrix_ned_to_cg(azimut,theta, phi):##
     return np.array([[ math.cos(azimut)* math.cos(theta),     math.sin(azimut)*math.cos(theta), -1.*math.sin(theta) ], 
                     [math.cos(azimut)* math.sin(theta)* math.sin(phi)- math.sin(azimut)* math.cos(phi),math.sin(azimut)* math.sin(theta)* math.sin(phi)+math.cos(azimut)* math.cos(phi),math.sin(phi)*math.cos(theta)], 
                     [math.cos(azimut)* math.sin(theta)* math.cos(phi)+ math.sin(azimut)* math.sin(phi),math.sin(azimut)* math.sin(theta)* math.cos(phi)-math.cos(azimut)* math.sin(phi),math.cos(phi)*math.cos(theta) ]])##
   
def get_matrix_a_to_cg(alpha, beta):
    return np.array([[ math.cos(alpha)*math.cos(beta), -1*math.cos(alpha)*math.sin(beta), -1.*math.sin(alpha) ], 
                     [                 math.sin(beta),                    math.cos(beta),                  0. ], 
                     [ math.sin(alpha)*math.cos(beta), -1*math.sin(alpha)*math.sin(beta),     math.cos(alpha) ]])


def get_angles_aer_cg(vel_cg):
    vel_abs = math.sqrt(vel_cg[0]**2 + vel_cg[1]**2 + vel_cg[2]**2)
    beta = math.asin(vel_cg[1] / vel_abs)
    alpha = -1*math.acos(vel_cg[0] / math.cos(beta) / vel_abs)
    return alpha, beta

def calc_mu(omega, dictFstate, dictRotor):
    vel_cg=np.array([dictFstate['vel_x_g'],dictFstate['vel_y_g'], dictFstate['vel_z_g']])
    vel_abs = math.sqrt(dictFstate['vel_x_g']**2 + dictFstate['vel_y_g']**2 + dictFstate['vel_z_g']**2)
    return vel_abs * math.cos(get_angles_aer_cg(vel_cg)[0]) / (omega * dictRotor['radius'])

def calc_mu_z(omega, dictFstate, dictRotor):
    vel_cg=np.array([dictFstate['vel_x_g'],dictFstate['vel_y_g'], dictFstate['vel_z_g']])
    vel_abs = math.sqrt(dictFstate['vel_x_g']**2 + dictFstate['vel_y_g']**2 + dictFstate['vel_z_g']**2)
    return -1. * vel_abs * math.sin(get_angles_aer_cg(vel_cg)[0]) / (omega * dictRotor['radius'])

def lambda_hover(omega, dictFstate, dictRotor):
    return vel_induced(omega, dictFstate, dictRotor) / (omega * dictRotor['radius'])

def inflow_equation(lamb_ind, omega, fstate, comp):
    mu = calc_mu(omega, fstate, comp)
    mu_z =calc_mu_z(omega, fstate, comp)
    lamb_hov = lambda_hover(omega, fstate, comp)
    return lamb_ind * math.sqrt( mu**2 + (mu_z + lamb_ind)**2 ) - lamb_hov**2

# =============================================================================
# def inflow_equation(lamb_ind, omega, fstate,alpha,dictRotor):##
#     v_ind=lamb_ind*omega*dictRotor['radius']##
#     mu =2*v_ind*math.cos(alpha)/(omega*dictRotor['radius']) ###
#     mu_z =-1*mu*math.tan(alpha) ###
#     ct=dictRotor['ct_sig']*dictRotor['n_bld']*dictRotor['chord']/(math.pi*dictRotor['radius'])##
#     lamb_hov =math.sqrt(ct/2)###
#     return lamb_ind * math.sqrt( mu**2 + (mu_z + lamb_ind)**2 ) - lamb_hov**2
# =============================================================================

def sum_forces_moments_cg(sol, fstate, components, results):
    
    forces = []
    moments = []
    powers_ind = []
    powers_prf = []
    powers_par = []
    
    #fuse = [comp for comp in components if comp['type'] == 'fuselage'][0]
    fstate['theta'] = sol[4]
    fstate['phi'] = sol[5]
    for comp in components:
        if comp['type'] == 'propeller':
            i = comp['index']

            # --- solve uniform inflow equation
            powerInd = power_induced(sol[i], fstate, comp)
            powerPrf = power_profile(sol[i], fstate, comp)###
            torque = np.array([ 0., 0., comp['rot'] * sol[i] / (powerInd + powerPrf) ])##
            # --- sum force and moment vectors and append to results list 
            powers_ind.append(powerInd)
            powers_prf.append(powerPrf)
            thrust = thrust_from_omega(sol[i], fstate, comp)###
            forces.append(np.array([0., 0., -1. * thrust ]))
            moments.append(np.cross(forces[-1], comp['pos']) + torque)
            #results['forces'] = sum(forces)##
            #results['moments'] = sum(moments)##
        elif comp['type'] == 'fuselage':
            gewicht_cg=np.matmul([0,0,calc_weight(comp['mass'])],get_matrix_ned_to_cg(0,sol[4],sol[5]))# --- set gravity vector and transform to cg system
            velocity_cg=np.matmul([fstate['vel_x_g'],fstate['vel_y_g'],fstate['vel_z_g']],get_matrix_g_to_cg(sol[4],sol[5]))# --- save velocity given in erdlotfestes system and transform to cg system
            [alpha,beta]=get_angles_aer_cg(velocity_cg)##alpha rechnen
            
            D= calc_drag_aer(fstate,comp)
            powerPar=-D* fstate['vel_x_g']
            powers_par.append(powerPar)
            drag_cg=np.matmul([-D,0,0], get_matrix_a_to_cg(alpha,beta))# --- calculate drag in aero system (e.g. use leverage according to concept) and transform to cg system
            forces_total=gewicht_cg+ drag_cg# --- sum force and moment vectors and append to results list 
            forces.append(np.array(forces_total))##
           # moments.append(np.array([0,0,0]))##
            
            M_N = np.array([0.,  D* 0.2 * alpha, 0.]) ### Hier sollte um ein passenderes Modell erweitert werden. 
            moments.append(M_N)
            # results['forces'] = sum(forces)##
            # results['moments'] = sum(moments)##

    residual = np.concatenate((sum(forces), sum(moments)), axis=None)
    
    results['power_ind'] = sum(powers_ind)
    results['power_prf'] = sum(powers_prf)
    results['power_par'] = sum(powers_par)
    
    return residual


#%%============================================================================
#    M A I N 
#%============================================================================= 
if __name__ == '__main__':
    
    # Define components of aircraft
    quadBody = {'type': 'fuselage',
                'mass': 1500.,
                'efpa': 1.25}   
    
    propFR = {'type':'propeller','index': 0,
              'radius': 2.0, 'n_bld': 3, 'ct_sig': 1.05, 'chord': 0.18, 'c_d0': 0.011,
              'pos': np.array([ 2.0, 2.0, 0.0 ]),'rot':1.}
    
    propFL = {'type':'propeller', 'index': 1,       
              'radius': 2.0, 'n_bld': 3, 'ct_sig': 1.05, 'chord': 0.18, 'c_d0': 0.011,
              'pos': np.array([ 2.0, -2.0, 0.0 ]),'rot':-1.}
    
    propAR = {'type':'propeller', 'index': 2,
              'radius': 2.0, 'n_bld': 3, 'ct_sig': 1.05, 'chord': 0.18, 'c_d0': 0.011,
              'pos': np.array([ -2.0, 2.0, 0.0 ]), 'rot':-1.}
    
    propAL = {'type':'propeller', 'index': 3,
              'radius': 2.0, 'n_bld': 3, 'ct_sig': 1.05, 'chord': 0.18, 'c_d0': 0.011,
              'pos': np.array([ -2.0, -2.0, 0.0 ]),'rot':1.}
    
    #ct = Schubbeiwert
    #sigma = Flächendichte 
    #ct_sigma = ct/sigma = Blattbelastung 
    #n_bld = Anzahl der Blätter 
    #rot = Rotationsrichtung 
    #Nehme für den Anfang an, dass dein Propeller ein Rotor ist
        
    # set x0 values
    omegaFR = 50.      # rad/s
    omegaFL = 50.      # rad/s
    omegaAR = 50.      # rad/s
    omegaAL = 50.      # rad/s
    theta_deg = 0.     # aircraft pitch angle (relative to erdlotfestes and thus ned cosy)
    phi_deg = 0.       # bank/roll angle (relative to erdlotfestes and thus ned cosy)
    theta = theta_deg * math.pi / 180.         # rad
    phi = phi_deg * math.pi / 180.             # rad
    
    # Initialize vector (list) for x0
    
    x0 = [omegaFR,omegaFL, omegaAR, omegaAL, theta, phi]
    
    components = [quadBody, propFR, propFL, propAR, propAL]
 
#%%============================================================================
#    C A S E: FLIGHT STATE 
#%=============================================================================

#    # Define flightstate
#    flightstate = {'vel_x_g': 50.,
#                   'vel_y_g': 0.,
#                   'vel_z_g': 0.,
#                   'dens': 1.225}
#
#    results = {} 
#    sol = optimize.root(sum_forces_moments_cg, x0, args=(flightstate, components, results))
#    print('\nSolution for hover found at:')
#    print('OmegaFR [rad/s]', sol.x[0])
#    print('OmegaFL [rad/s]', sol.x[1])
#    print('OmegaAR [rad/s]', sol.x[2])
#    print('OmegaAL [rad/s]', sol.x[3])
#    print('Theta [deg]', sol.x[4]*180/math.pi)
#    print('Phi [deg]', sol.x[5]*180/math.pi)
#    sys.exit()
    
    
#%%============================================================================
#    C A S E: FLIGHT ENVELOPE 
#%=============================================================================    
    
    # Define flightstate
    flightstate = {'vel_x_g': 0.,
                   'vel_y_g': 0.,
                   'vel_z_g': 0.,
                   'dens': 1.225}
    
    #sum_forces_moments_cg(x0, dens, components)
    velocities = range(1,80)
    #velocities = range(15,16)
    
    omegasA = []
    omegasF = []
    thetas = []
    pInd = []
    pPrf = []
    pPar = []
    for vel in velocities:
        flightstate.update({'vel_x_g': vel})
        results = {}
        sol = optimize.root(sum_forces_moments_cg, x0, args=(flightstate, components, results))
        omegasF.append(sol.x[0])
        omegasA.append(sol.x[2])
        thetas.append(sol.x[4]*180/math.pi)
        pInd.append(results['power_ind']/1000.)
        pPrf.append(results['power_prf']/1000.)
        pPar.append(results['power_par']/-1000.)
        pTot=np.array(pInd)+np.array(pPrf)+np.array(pPar)
        

    fig0, ax0 = plt.subplots()
    ax0.plot(velocities, thetas, label='Theta')
    ax0.legend()
    ax0.grid()
    plt.xlabel('Velocity v_x')
    plt.ylabel('Flight Angle Theta')
    plt.show()
    
    fig1, ax1 = plt.subplots()
    ax1.plot(velocities, omegasF, label='Front')
    ax1.plot(velocities, omegasA, label='Back')
    plt.xlabel('Velocity v_x')
    plt.ylabel('RPM')
    plt.legend()
    plt.grid()
    plt.show()
    
    fig2, ax2 = plt.subplots()
    ax2.plot(velocities, pInd, label='pInd')
    ax2.plot(velocities, pPrf, label='pPrf')
    ax2.plot(velocities, pPar, label='pPar')
    ax2.plot(velocities, pTot, label='pTot')
    plt.xlabel('Velocity v_x')
    plt.ylabel('Power')
    plt.legend()
    plt.grid()
    plt.show()
        
    sys.exit()    


