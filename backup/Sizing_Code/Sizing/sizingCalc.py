#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 13:46:10 2020

@author: christian
"""

from scipy import optimize
from calculation_Hover import getPowerHover
from calculation_Cruise import getPowerCruise
import numpy as np

def funcComponents(m, dData):
    
    dData['masses'] = {}
    dData['results'] = {}
  
    vCrs = dData['mission']['cruiseVelo']             # m/s
    vClimb = dData['mission']['climbVelo']              # m/s

    tRes = dData['mission']['timeReserve']            # m
    rCrs = dData['mission']['cruiseRange']            # m
    hAGL = dData['mission']['altitudeAGL']
    
    # Calc Powers
    dData['results']['powers'] = {}
    dData = getPowerCruise(m, dData)
    pCrs = dData['results']['powers']['cruisePower']


    dData = getPowerHover(m, dData)
    pClimb = dData['results']['powers']['climbPower']
    pSink = dData['results']['powers']['sinkPower']
    
    tCrs = rCrs / vCrs
    tClimb = hAGL / vClimb
    tSink = hAGL / (vClimb*1.5)
    
    eta = dData['battery']['eta']                        # -
    
    # Hydrogen Fuel Cell:
    specPowerFC = dData['h2']['specPowerFC']
    specEnergyH2 = dData['h2']['specEnergyH2']
    etaFC = dData['h2']['etaFC']
    H2Proportion = dData['h2']['H2Proportion']
    
    fP = 1.1                        # factor power FC
    
    FCPower = pCrs * fP / eta
    
    fuelCell = FCPower / specPowerFC
    fuel = FCPower * ((tCrs*2 + tRes) / fP + (tClimb + tSink)*2) / (specEnergyH2 * etaFC)
    tank = fuel / H2Proportion - fuel
    
    dData['results']['fuelNeeded'] = FCPower * (tCrs*2 / fP + (tClimb + tSink)*2) / (specEnergyH2 * etaFC)

    # Battery:
    voltage = dData['battery']['voltage']                # V                                 # s
    specEnergyBat = dData['battery']['specEnergy'] * 3600.      # Ws/kg
    cMax = dData['battery']['cRate'] / 3600.          # 1/s
    motorError = dData['battery']['motorError']        
    # ---------------------------------------------------------------------
    reqCap = capacityBat(pClimb-FCPower*eta, pSink-FCPower*eta, tClimb, tSink, voltage, eta, cMax)
    reqMassEnergy = reqCap * voltage / specEnergyBat
    # ---------------------------------------------------------------------
    dData['results']['batteryEnergy'] = reqCap * voltage
    dData['results']['batteryPower'] = reqCap * voltage * cMax
    
    dData['results']['FCPower'] = FCPower
    
    energyClimb = (pClimb-FCPower*eta) / eta / voltage * tClimb
    energySink = (pSink-FCPower*eta) / eta / voltage * tSink
    dData['results']['energyNeeded'] = (energyClimb + energySink) * 2
    
    battery = reqMassEnergy
    ptrain = motorModelNg(max(pClimb, pCrs)) * motorError
        
    mPayload = dData['mission']['nPax'] * 88
    mCrew = dData['mission']['massCrew']
    fStruct = dData['design']['fStruct'] 
    fOth = dData['design']['fOth']
    useful = mPayload + mCrew
    struct = fStruct * m[0]
    other = useful * fOth
    
    dData['masses']['struct'] = struct
    dData['masses']['ptrain'] = ptrain
    dData['masses']['fuel'] = fuel    
    dData['masses']['useful'] = useful
    dData['masses']['other'] = other
    dData['masses']['dgm'] = m
    dData['masses']['tank'] = tank
    dData['masses']['battery'] = battery
    dData['masses']['fuelCell'] = fuelCell
    
    return dData

def motorModelNg(pReq): # based on state of the art electric motors for aviation
    return np.exp(-0.89 + 0.89 * np.log(pReq / 1000.))

def capacityBat(powerClimb, powerSink, tClimb, tSink, voltage, eta, cMax):
    energyClimb = powerClimb / eta / voltage * tClimb
    energySink = powerSink / eta / voltage * tSink
    energywise = (energyClimb + energySink) / 0.8        # safety factor of 1/0.8
    powerwise = powerClimb / eta / voltage / cMax * 1.3  # safety faktor of 1.3
    capEnergy = max(energywise, powerwise)
    return capEnergy

def objFuncSolveMass(m, dArgs):
    
    dArgs = funcComponents(m, dArgs)    
    struct = dArgs['masses']['struct']
    ptrain = dArgs['masses']['ptrain']
    fuel = dArgs['masses']['fuel']
    fuelCell = dArgs['masses']['fuelCell']
    tank = dArgs['masses']['tank']
    battery = dArgs['masses']['battery']
    useful = dArgs['masses']['useful']
    other = dArgs['masses']['other']
    
    residual =  m - struct - ptrain - fuel - battery - fuelCell - tank - useful - other
    return residual

def getAircraftMass(dData):
    res = 0.
    startDGW = (dData['mission']['massCrew'] + dData['mission']['nPax']*88) / 0.25    

    res = optimize.root(objFuncSolveMass, startDGW, args=dData)
    if res.success:
        dData['masses']['dgm'] = res.x[0]
        dData['info'] = True
    else:
        dData['info'] = False
    
    return dData

def exactMass(dData):
    return getAircraftMass(dData)['masses']['dgm']