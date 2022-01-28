#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 15:54:07 2021

@author: sumeetkumar
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
import numpy as np
##############################################################################
############## Quantities corresponding to a sample helicopter 
##############################################################################
Vinf_max = 100    #max. forward speed
k = 1.15    #kappa
W = 3400*9.8
d_area = 1.1   #drag area
R = 5.5
Vtip = 210
sigma = 0.08
Cd0 = 0.01

p0 = 101325 #std atmospheric pressure
T0 = 288    #std atmposheric temp
R_gas = 287.04  #gas constant
Ae = 40
Be = 0.25
h = np.array([-2000,0,2000,4000,6000])
N_eng = 1   #no of engines
##############################################################################
################### derived quantities
##############################################################################
A = np.pi*R**2
Vinf = np.linspace(0,Vinf_max,41)
mu = Vinf/Vtip
##############################################################################
################## functions to calculate relevant quantities  
##############################################################################
def get_vi(Vinf,T=W,alt=0):    
    
    vh = (T/(2*get_rho(alt=alt)*A))**0.5    
    vi_arr = (-0.5*Vinf**2 + (vh**4 + 0.25*Vinf**4)**0.5)**0.5
        
    return vi_arr

def get_Pi_id(Vinf,T=W,alt=0):
    
    return T*get_vi(Vinf,alt=alt)

def get_P0(Vinf,alt=0):
    
    return 0.125*get_rho(alt=alt)*A*(Vtip**3)*sigma*Cd0*(1+4.65*((Vinf/Vtip)**2))
 
def get_Pp(Vinf,alt=0):
    
    return 0.5*get_rho(alt=alt)*(Vinf**3)*d_area

def get_P(Vinf,alt=0):   #default result at 0 altitude

    Pi = k*get_Pi_id(Vinf, alt=alt)
    P0 = get_P0(Vinf, alt=alt)
    Pp = get_Pp(Vinf, alt=alt) 
    
    return Pi+P0+Pp

def get_rho(alt=0):
    
    return get_p(alt)/(R_gas*get_T(alt))

def get_p(alt):
    
    return p0*(1-0.0065*alt/T0)**5.2561

def get_T(alt):
    
    return T0*(1-0.0065*alt/T0)
##############################################################################
##############    Plotting and verifying the power curves     ################
##############################################################################
Pi = k*get_Pi_id(Vinf)
P0 = get_P0(Vinf)
Pp = get_Pp(Vinf)  
P = get_P(Vinf)

plt.figure()
plt.plot(mu, P/1000, color="tomato", linestyle="-", label="P")
plt.plot(mu, Pi/1000, color="olivedrab", linestyle="-", label="Pi")
plt.plot(mu, P0/1000, color="firebrick", linestyle="-", label="P0")
plt.plot(mu, Pp/1000, color="cornflowerblue", linestyle="-", label="Pp")

plt.ylim((0, 900))
plt.xlim((0, 0.5))
plt.ylabel(r"Power [kW]")
plt.xlabel(r"$\mu$")
plt.legend()
plt.savefig('Power.png')
##############################################################################
##################     Plotting delta and theta curves     ###################
##############################################################################
p = get_p(h)
T = get_T(h)
delta = p/p0
theta = T/T0

plt.figure()
plt.plot(h, delta, color="red", linestyle="-", label=r"$\delta$")
plt.plot(h, theta, color="blue", linestyle="-", label=r"$\theta$")
plt.plot(h, theta**0.5, color="green", linestyle="-", label=r"$\sqrt{\theta}$")
plt.plot(h, delta*theta**0.5, color="black", linestyle="-", label=r"$\delta\sqrt{\theta}$")

plt.xlabel("h [m]")
plt.legend()
plt.savefig('delta_theta.png')
##############################################################################
####  Plotting fuel consumption and power curves for different altitudes  ####
##############################################################################

color_scheme1 = ["olivedrab","firebrick","cornflowerblue","lightgray","darkgrey"] #same size as h
color_scheme2 = ["red","lime","deepskyblue","black","grey"] #same size as h

#Additional functions for creating W_fuel and Power vs altitude plots
def def_fig_axes():
    fig = plt.figure(figsize=(8,3.5))
    gs = gridspec.GridSpec(1, 2, figure=fig)
    ax_wfuel = fig.add_subplot(gs[0])
    ax_p_alt = fig.add_subplot(gs[1])
    
    return fig, ax_wfuel, ax_p_alt

def plot_wf_p(n_eng,ax_wfuel,ax_p_alt):
    Wfuel_dot = np.zeros((len(Vinf),len(h)))
    P_alt = np.zeros((len(Vinf),len(h)))
    for i,(alt,c) in enumerate(zip(h,color_scheme2)):    
        P_alt[:,i] = get_P(Vinf,alt=alt)
        Wfuel_dot[:,i] = n_eng*(Ae*delta[i]*(theta[i]**0.5) + Be*(P_alt[:,i]/1000)*(1/n_eng))
    
        ax_wfuel.plot(mu, Wfuel_dot[:,i], color=c)
        if n_eng==N_eng:   #to avoid repeating labels in legend
            labl=f"h={alt} m"
        else:
            labl=None                
        ax_p_alt.plot(mu, P_alt[:,i]/1000, color=c, label=labl)
    # ax_wfuel.set_title(f"Number of engines = 2")
  
def fig_details(ax_wfuel,ax_p_alt):
    ax_wfuel.set_ylabel(r"$\dot{W}_{fuel}$ [kg/h]")
    ax_wfuel.set_xlabel(r"$\mu$")
    ax_p_alt.set_ylabel(r"Power [kW]")
    ax_p_alt.set_xlabel(r"$\mu$")
    ax_wfuel.set_ylim(100,450)
    ax_p_alt.set_ylim(0,900)
    ax_p_alt.legend()

# Plotting and verifying the W_fuel and Power vs altitude curves
Fig, Ax_wfuel, Ax_p_alt = def_fig_axes()
plot_wf_p(N_eng,Ax_wfuel,Ax_p_alt)  
fig_details(Ax_wfuel, Ax_p_alt)

##############################################################################
#Creating a gif to see the effect of 'N_eng' on the 'W_fuel' curves
# from celluloid import Camera
# Fig, Ax_wfuel, Ax_p_alt = def_fig_axes()
# camera = Camera(Fig)

# for n_eng in range(1,N_eng+1,1):
#     if N_eng==1:
#         raise Exception(f"N_eng is currently {N_eng}. It should be greater than 1 to create a gif")
#     plot_wf_p(n_eng,Ax_wfuel,Ax_p_alt)
#     camera.snap()

# fig_details(Ax_wfuel, Ax_p_alt)
    
# animation = camera.animate()  
# animation.save('Wfuel.gif', writer = 'pillow')