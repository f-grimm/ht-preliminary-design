#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 13:46:10 2020

@author: christian
"""

from calculation_Hover import getPowerHover
from calculation_Cruise import getPowerCruise

def funcComponents(dData):
    
    struct = dData['masses']['struct']
    ptrain = dData['masses']['ptrain']
    fuel = dData['masses']['fuel']
    other = dData['masses']['other']
    tank = dData['masses']['tank']
    battery = dData['masses']['battery']
    fuelCell = dData['masses']['fuelCell']
    useful = dData['masses']['useful']
    
    m = useful + other + battery + fuelCell + fuel + tank + struct + ptrain
    
    vCrs = dData['mission']['cruiseVelo']             # m/s
    vClimb = dData['mission']['climbVelo']              # m/s
    
    tRes = dData['mission']['timeReserve']            # m
    hAGL = dData['mission']['altitudeAGL']
    
    dData = getPowerCruise(m, dData)
    pCrs = dData['results']['powers']['cruisePower']
    
    dData = getPowerHover(m, dData)
    pClimb = dData['results']['powers']['climbPower']
    pSink = dData['results']['powers']['sinkPower']
    
    tClimb = hAGL / vClimb
    tSink = hAGL / (vClimb*1.5)
    
    eta = dData['battery']['eta']                        # -
    
    # Battery:
    specEnergyBat = dData['battery']['specEnergy'] * 3600.      # Ws/kg
    cMax = dData['battery']['cRate'] / 3600.
    
    specEnergyH2 = dData['h2']['specEnergyH2']
    etaFC = dData['h2']['etaFC']
    H2E = fuel * specEnergyH2 * etaFC
    
    batE = battery * specEnergyBat * 0.8
    batP = battery * specEnergyBat * cMax / 1.3
    
    batP = min(batP, batE / (tClimb + tSink))

    energyH2Cruise = (H2E - ((pClimb / eta - batP) * tClimb + (pSink / eta - batP) * tSink) * 2)
        
    dData['results']['distance'] = vCrs * (energyH2Cruise * eta / pCrs - tRes) /2
    dData['results']['distance_double'] = vCrs * (energyH2Cruise * eta / pCrs - tRes)
    
    return dData