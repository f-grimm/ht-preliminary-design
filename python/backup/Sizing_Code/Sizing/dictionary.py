#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 13:41:12 2020

@author: christian
"""
def getDict(state):
    
    nCrew = 1
    minutesReserve = 20.
    mCrew = 92 * nCrew                # repr. for men
    tReserve = 60 * minutesReserve    # regulation: 20mins in fwd flight
    
    
    dRequiredMission = {'ident':          'Sizing',     
                        'cruiseRange':      40000.,     # m
                        'cruiseVelo':          55.,     # m/s - 250km/h (Web 300km/h)
                        'climbVelo':            5.,     # m/s
                        'massCrew':          mCrew,     # kg
                        'nPax':                  4,     # -
                        'altitude':          1000.,     # m ü NN
                        'altitudeAGL':        500.,
                        'timeReserve':    tReserve}     # s
    
    dDesignAircraft = {'ident':  'tilt-prop',
                       'fStruct':       0.25,         # - (mStruc = fStruc * DGW)
                       'fOth':           1.}         # - (mOth = mUse * fOth)
    
    if state == 'SoA':
        dHydrogen = {    'fueltype':             'hydrogen',
                         'eta':                       0.85,      # - (motor efficiency)
                         'specPowerFC':           1.5*1000,      # W/kg ()
                         'specEnergyH2':   33.33*1000*3600,      # J/kg (specific energy H2)
                         'etaFC':                      0.4,      # - (Fuel Cell efficiency)
                         'H2Proportion':             0.055       # - (Proportion of H2 of the H2-Tank-System)
                        }               
        
        dBatteryLiIon = {'fueltype':      'battery',
                         'nCycles':            1000,            # -
                         'specEnergy':      108/1.2,            # Wh/kg (pack specific = cell specific / 1.2) 
                         'eta':                0.85,            # - (motor efficiency)
                         'voltage':             2.2,            # V 425 - 650 but doesnt have a direct influence on weight if C-rate arbitrary
                         'dcp':                1.03,            # - (battery discharge parameter - see L.Traub )
                         'motorError':            1,            # - (factor multiplied on the electric motor weight)
                         'cRate':               20.,            # - (max discharge rate)
                         'fDoD':                0.8}
    
    elif state == 'future':
        dHydrogen = {    'fueltype':            'hydrogen',
                         'eta':                       0.85,      # - (motor efficiency)
                         'specPowerFC':            2.*1000,       # W/kg ()
                         'specEnergyH2':   33.33*1000*3600,      # J/kg (specific energy H2)
                         'etaFC':                     0.45,      # - (Fuel Cell efficiency)
                         'H2Proportion':             0.065       # - (Proportion of H2 of the H2-Tank-System)
                         }
        
        dBatteryLiIon = {'fueltype':      'battery',
                         'nCycles':            1000,             # -
                         'specEnergy':        140/1.2,           # Wh/kg (pack specific = cell specific / 1.2) 
                         'eta':                0.85,             # - (motor efficiency)
                         'voltage':             2.2,             # V 425 - 650 but doesnt have a direct influence on weight if C-rate arbitrary
                         'dcp':                1.03,             # - (battery discharge parameter - see L.Traub )
                         'motorError':            1,             # - (factor multiplied on the electric motor weight)
                         'cRate':               20.,             # - (max discharge rate)
                         'fDoD':                0.8}
       
    ### Configuration
       
    body = {'type': 'fuselage', 'efpa': 1.25}
    tiltProp = {'type':'propeller',
                'radius': 0.6, 'n_bld': 5, 'ct_sig': 0.15, 'chord': 0.18, 'c_d0': 0.011}
    hoverProp = {'type':'propeller',
                'radius': 0.6, 'n_bld': 5, 'ct_sig': 0.15, 'chord': 0.18, 'c_d0': 0.011}
    
    # Hover Configuration:
    props = []
    for n in range(0,4):
        props.append(hoverProp)
    for n in range(0,12):
        props.append(tiltProp)
    hover = [body] + props
    
    # Cruise Configuration:
    props = []
    wingF = {'type':'wing',
                 'chord': 0.8, 'l': 6., 'e': 0.85}
    
    wingB = {'type':'wing',
                'chord': 1., 'l': 8., 'e': 0.85}
    
    for n in range(0,12):
        props.append(tiltProp)
    cruise = [body, wingF, wingB] + props
    
    dConfiguration = {'hover': hover, 'cruise': cruise}
        
    dSizingCase = {'design': dDesignAircraft, 'mission': dRequiredMission,
                   'h2': dHydrogen, 'battery': dBatteryLiIon, 'configuration': dConfiguration}
    
    return dSizingCase