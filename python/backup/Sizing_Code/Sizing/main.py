#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 14:20:19 2019

Nicolas André, MSc
Institute of Helicopter Technology (HT)
Techical University of Munich (TUM)

NOTICE:  The code provided herein is solely for the purpose of the 
project seminar "Urban Air Mobility" at the Institute of Helicopter Technology.
All information contained herein is, and remains
the property of the Institute of Helicopter Technology.
The intellectual and technical concepts contained
herein are proprietary to Institute of Helicopter Technology, 
and may be protected by copyright law.
Dissemination of this information or reproduction of this material
is strictly forbidden unless prior written permission is obtained
from Institute of Helicopter Technology.

@author: gu95pug
"""

# used for sizing
from sizingCalc import getAircraftMass, exactMass

# used for analysis and plotting
from analysis import calcSensitivitiesOAT

# used for carpet and plotting
import matplotlib.pyplot as plt
import numpy as np

# used to get the configuration dict
from dictionary import getDict

# used for reverse range calc
import rangeCalc

# used to calc powers at constant dgm
from calculation_Cruise import stationaryCrs
from calculation_Hover import stationaryHov

# %%============================================================================
#     S I Z I N G
# %=============================================================================
# ---- sizing
dSizingCase = getDict('SoA')
dResult = getAircraftMass(dSizingCase) # (deterministic)
if dResult['info'] == True:
    dgmAircraft = dResult['masses']['dgm']
    print('\n############# Sizing Baseline #################')      
    print('The Aircraft Design Gross Weight at given \nRequirements is DGW = {:6.1f} kg'.format(dgmAircraft))
    
    print('\n############# All masses #################')  
    for k, v in dResult['masses'].items():
        if k != 'dgm':
            print(' {:14s} {:10.1f} kg'.format(k, v))
        
    print('\n############# Power Requirements #################')    
    for k, v in dResult['results']['powers'].items():
        print(' {:14s} {:10.0f} kW'.format(k, v/1000.))
        
    print('\n############# Mission Details #################')   
    for k, v in dResult['mission'].items():
        print(' {:16s} {}'.format(k, v)) 
        
    print('\n############# Results #################')  
    
    batEn = dResult['results']['batteryEnergy']/1000./3600.
    batP = dResult['results']['batteryPower']/1000
    fcP = dResult['results']['FCPower']/1000.
    print(' {:16s} {:.0f} kWh'.format("Battery Energy", batEn))
    print(' {:16s} {:.0f} kW'.format("Battery Power", batP))
    print(' {:16s} {:.0f} kW'.format("Fuel Cell Power", fcP))
    print(' {:16s} {:.0f} kg'.format("H2 needed", dResult['results']['fuelNeeded']))
    print(' {:16s} {:.0f} kWh'.format("BatEnergy needed", dResult['results']['energyNeeded']/1000./3600.))
    
else:
    print("No solution")
    
# %%============================================================================
#     S T A T I O N A R Y  P O W E R S
# %=============================================================================
stationaryCrs(dSizingCase)
stationaryHov(dSizingCase)

# %%============================================================================
#     M A S S E S  P I E  C H A R T
# %=============================================================================
xDict = dSizingCase['masses'].copy()

xDict.pop('dgm', None)
values = list(xDict.values())
labels = list(xDict.keys())


labels = [labels for _,labels in sorted(zip(values,labels))]
labels = ['power train' if x=='ptrain' else x for x in labels]
labels = ['fuel cell' if x=='fuelCell' else x for x in labels]
values.sort()
fig1, ax1 = plt.subplots()

ax1.pie(values, labels=labels, autopct='%1.1f%%', startangle=180)
ax1.axis('equal')
plt.title('Mass Distribution')
plt.show()

# %%============================================================================
#     S I Z I N G ( V E L O C I T I E S )
# %=============================================================================

velocities = np.arange(30,71, 1.)

masses = []
c_ls = []
ltods = []
masses_f = []
c_ls_f = []
ltods_f = []
crsPowers = []
hovPowers = []

for vel in velocities:
    dData = dSizingCase.copy()
    dData['mission']['cruiseVelo'] = vel
    dResult = getAircraftMass(dData)
    if dResult['info'] == True:
        masses.append(dResult['masses']['dgm'])
        c_ls.append(dResult['results']['c_l'])
        ltods.append(dResult['results']['ltod'])
        crsPowers.append(dResult['results']['powers']['cruisePower'] / 1000.)
        hovPowers.append(dResult['results']['powers']['climbPower'] / 1000.)
    else:
        masses.append(None)
        c_ls.append(None)
        ltods.append(None)
        crsPowers.append(None)
        hovPowers.append(None)
        
    dData = getDict('future')
    dData['mission']['cruiseVelo'] = vel
    dResult = getAircraftMass(dData)
    if dResult['info'] == True:
        masses_f.append(dResult['masses']['dgm'])
        c_ls_f.append(dResult['results']['c_l'])
        ltods_f.append(dResult['results']['ltod'])
    else:
        masses_f.append(None)
        c_ls_f.append(None)
        ltods_f.append(None)

plt.title('Sizing over velocity (range: 40km)')
plt.plot(velocities, masses, label='SoA')
plt.plot(velocities, masses_f, label='future')
plt.plot(55, dgmAircraft, 'o', color='red')
plt.axis([30, 70, 1500, 4100])
plt.legend()
plt.xlabel('Velocity (m/s)')
plt.ylabel('Mass (kg)')
plt.grid()
plt.show()

plt.title('Sizing over velocity (range: 40km)')
plt.plot(velocities, c_ls, label='SoA')
plt.plot(velocities, c_ls_f, label='future')
plt.axhline(y=1.6, color='black', linestyle="--")
plt.axis([30, 70, 0, 3])
plt.legend()
plt.xlabel('Velocity (m/s)')
plt.ylabel('C_l')
plt.grid()
plt.show()

plt.title('Sizing over velocity (range: 40km)')
plt.plot(velocities, ltods, label='SoA')
plt.plot(velocities, ltods_f, label='future')
plt.axis([30, 70, 0, 9])
plt.legend()
plt.xlabel('Velocity (m/s)')
plt.ylabel('L/D')
plt.grid()
plt.show()

plt.title('Sizing over velocity (range: 40km)')
plt.plot(velocities, hovPowers, label='Climb')
plt.plot(velocities, crsPowers, label='Cruise')
plt.axis([30, 70, 0, 1000])
plt.legend()
plt.xlabel('Velocity (m/s)')
plt.ylabel('Power (kW)')
plt.grid()
plt.show()
    
# %%============================================================================
#     S I Z I N G ( R A N G E )
# %=============================================================================

distances = np.arange(25,215, 5.)

masses = []
masses_f = []


for dist in distances:
    dData = getDict('SoA')
    dData['mission']['cruiseRange'] = dist * 1000
    dData = getAircraftMass(dData)
    if dData['info'] == True:
        masses.append(dData['masses']['dgm'])
    else:
        masses.append(None)
        
    dData = getDict('future')
    dData['mission']['cruiseRange'] = dist * 1000
    dData = getAircraftMass(dData)
    if dData['info'] == True:
        masses_f.append(dData['masses']['dgm'])
    else:
        masses_f.append(None)
        
plt.title('Sizing over range (V: 55 m/s)')
plt.axhline(3175, color='black', linestyle='--')
plt.plot(distances, masses, label='SoA')
plt.plot(distances, masses_f, label='future')
plt.plot(40, dgmAircraft, 'o', color='red')
plt.axis([25, 210, 1500, 4100])
plt.legend()
plt.xlabel('Range (km)')
plt.ylabel('Mass (kg)')
plt.grid()
plt.show()
    
    
# %%============================================================================
#     R A N G E
# %=============================================================================

usefulmasses = np.arange(0,510, 10.)
distances = []
distances_double = []

for useful in usefulmasses:
    
    dData = dResult.copy()
    
    dData['masses']['useful'] = useful
    dData = rangeCalc.funcComponents(dData)
    
    distances.append(dData['results']['distance'] / 1000.)
    distances_double.append(dData['results']['distance_double'] / 1000.)
    

plt.title('Variation of the useful mass (sized for 444kg)')
plt.plot(usefulmasses, distances, label='with return')
plt.plot(usefulmasses, distances_double, label='no return')
# plt.axis([0, 500, 35, 120])
plt.legend()
plt.xlabel('Useful Mass (kg)')
plt.ylabel('Range (km)')
plt.grid()
plt.show()

# %%============================================================================
#     S E N S I T I V I T I E S
# %=============================================================================
    
dSizingCase = getDict('SoA')
iSamples = 25

lower = -0.12
bLow = 1.0+lower
upper = 0.12
bUp = 1.0+upper

interval = np.linspace(lower,upper,iSamples)

dVariation = {'battery':
                  {'specEnergy': np.linspace(dSizingCase['battery']['specEnergy'] * bLow, dSizingCase['battery']['specEnergy'] * bUp, iSamples)},
              'h2':
                  {'etaFC': np.linspace(dSizingCase['h2']['etaFC'] * bLow, dSizingCase['h2']['etaFC'] * bUp, iSamples),
                    'specPowerFC': np.linspace(dSizingCase['h2']['specPowerFC'] * bLow, dSizingCase['h2']['specPowerFC'] * bUp, iSamples)
                  },
            'design':
                {'fStruct': np.linspace(dSizingCase['design']['fStruct'] * bLow, dSizingCase['design']['fStruct'] * bUp, iSamples),
                  'fOth': np.linspace(dSizingCase['design']['fOth'] * bLow, dSizingCase['design']['fOth'] * bUp, iSamples)
                    }}

dSensAnalysis = calcSensitivitiesOAT(exactMass, dSizingCase, dVariation)

figSensAnalysis = plt.figure(figsize=(6,4))

plt.title('One-at-a-time Sensitivity Analysis', fontsize=12)

for k0, v0 in dSensAnalysis.items():
    for k1, v1, in v0.items():
        plt.plot(interval, v1, label=k1, linewidth=0.8)

plt.xlabel('Variation from baseline')
plt.ylabel('Mass (kg)')
plt.grid()
labels = ['Spec. Energy Battery', 'Fuel Cell Efficiency', 'Spec. Power Fuel Cell', 'Struct Ratio', 'Other Masses Ratio']
plt.legend(labels)
plt.show()

# %%============================================================================
#     S I Z I N G ( F U T U R E )
# %=============================================================================

velocities = np.arange(30,70, 1.)

masses = []

for vel in velocities:
    dData = getDict('future')
    dData['mission']['cruiseVelo'] = vel
    dData['mission']['cruiseRange'] = 150 * 1000.
    dResult = getAircraftMass(dData)
    if dResult['info'] == True:
        masses.append(dResult['masses']['dgm'])
        if vel == 55.:
            print(' Cruise Power with 150km range at 55m/s: {:6.0f}kW'.format(dResult['results']['powers']['cruisePower']/1000.))
            print(' Climb Power with 150km range at 55m/s: {:6.0f}kW'.format(dResult['results']['powers']['climbPower']/1000.))
    else:
        masses.append(None)

plt.title('Future with Range of {}km'.format(dData['mission']['cruiseRange']/1000.))
plt.plot(velocities, masses)
plt.axhline(3175, color='black', linestyle='--')
plt.axis([30, 70, 2700, 3500])
plt.xlabel('Velocity (m/s)')
plt.ylabel('Mass (kg)')
plt.grid()
plt.show()

masses = []

for vel in velocities:
    dData = getDict('future')
    dData['mission']['cruiseVelo'] = vel
    dData['mission']['nPax'] = 6
    dResult = getAircraftMass(dData)
    if dResult['info'] == True:
        masses.append(dResult['masses']['dgm'])
        if vel == 55.:
            print(' Cruise Power with {:.0f} Pax at 55m/s: {:6.0f}kW'.format(dData['mission']['nPax'],dResult['results']['powers']['cruisePower']/1000.))
            print(' Climb Power with {:.0f} Pax at 55m/s: {:6.0f}kW'.format(dData['mission']['nPax'],dResult['results']['powers']['climbPower']/1000.))
    else:
        masses.append(None)

plt.title('Future with {} Passengers'.format(dData['mission']['nPax']))
plt.plot(velocities, masses)
plt.axhline(3175, color='black', linestyle='--')
plt.axis([30, 70, 2700, 3500])
plt.xlabel('Velocity (m/s)')
plt.ylabel('Mass (kg)')
plt.grid()
plt.show()

