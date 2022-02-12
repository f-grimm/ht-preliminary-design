#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 11:38:48 2020

@author: christian
"""

import matplotlib.pyplot as plt

def generatePlot(fstate, var, values, dData, results):
    
    dgm = dData['masses']['dgm']
    
    powers = results['powers']
    pInd = powers['power_ind']
    pPrf = powers['power_prf']
    pGes = powers['power_ges']

    if fstate['type'] == 'cruise':
        
        pMax = dData['results']['FCPower'] * dData['h2']['eta'] / 1000.
        
        forces = results['forces']
        
        pPar = powers['power_par']
        c_ls = results['c_ls']
        drags_ind = forces['drags_ind']
        drags_par = forces['drags_par']
        drags_ges = forces['drags_ges']
        ltodrags = forces['ltodrags']
        
        plt.title('Lift Coefficent (dgm: {:.0f}kg)'.format(dgm), fontsize=12)
        plt.plot(values, c_ls)
        plt.axis([30, 70, 0, 3])
        plt.axhline(y=1.6, color='black', linestyle="--")
        plt.xlabel('Velocity (m/s)')
        plt.ylabel('C_l')
        plt.grid()
        plt.show()
        
        plt.title('Drag Forces (dgm: {:.0f}kg)'.format(dgm), fontsize=12)
        plt.plot(values, drags_ind, label='Ind')
        plt.plot(values, drags_par, label='Par')
        plt.plot(values, drags_ges, label='Ges')
        plt.axis([30, 70, 0, 5])
        plt.legend()
        plt.xlabel('Velocity (m/s)')
        plt.ylabel('Drag (kN)')
        plt.grid()
        plt.show()
        
        plt.title('Cruise Power (dgm: {:.0f}kg)'.format(dgm), fontsize=12)
        plt.plot(values, pInd, label='Induced')
        plt.plot(values, pPrf, label='Profile')
        plt.plot(values, pPar, label='Parasitic')
        plt.plot(values, pGes, label='Total')
        plt.axhline(y=pMax, color='black', linestyle="--")
        plt.xlabel('Velocity (m/s)')
        plt.ylabel('Power (kW)')
        plt.axis([30, 70, 0, 350])
        plt.legend()
        plt.grid()
        plt.show()
        
        plt.title('Lift to Drag Ratio (dgm: {:.0f}kg)'.format(dgm), fontsize=12)
        plt.plot(values, ltodrags)
        plt.xlabel('Velocity (m/s)')
        plt.ylabel('L/D')
        plt.grid()
        plt.show()
    
    elif fstate['type'] == 'hover':
        
        pMax = (dData['results']['batteryPower'] + dData['results']['FCPower']) * dData['battery']['eta'] / 1000.
        
        BatFailP = (dData['results']['batteryPower']/2 + dData['results']['FCPower']) * dData['battery']['eta'] / 1000.
        
        pClb = powers['power_clb']
        
        plt.title('Climb/Sink (dgm: {:.0f}kg)'.format(dgm), fontsize=12)
        plt.axvline(x=0, color='grey')
        plt.plot(values, pInd, label='Induced', color='blue')
        plt.plot(values, pPrf, label='Profile', color='orange')
        plt.plot(values, pClb, label='Climb', color='grey')
        plt.plot(values, pGes, label='Total', color ='brown')
        plt.axhline(y=pMax, color='black', linestyle="--")
        plt.axhline(y=BatFailP, color='red', linestyle="--")
        plt.axis([-10,10, -100, 900])
        plt.legend()
        plt.xlabel('Velocity (m/s)')
        plt.ylabel('Power (kW)')
        plt.grid()
        plt.show()