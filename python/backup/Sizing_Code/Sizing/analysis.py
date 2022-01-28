#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 17:53:09 2019

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
#test Kommentar
import numpy as np

def calcCarpet(fun, dCase, dCarpet):
    """Calculate a data for a two-dimensional carpet plot
    
    Parameters
    ----------
    fun : function
        The function that needs to be solved for Y
        Has to be of the form fun(dictionary) with dData
    dCase : dictionary
        What it does
    dCarpet : dictionary
        What it does

    Returns
    -------
    dol : list
        a list of strings used that are the header columns
    """
 
    print('\nInitiate carpet calculation...')
    keys0 = []
    keys1 = []
    for k0, v0 in dCarpet.items():
        for k1, v1 in v0.items():
            dim = len(v1)
            keys0.append(k0)
            keys1.append(k1)
            #print(k0, k1, v1)   
      
    # get the matrix X for x placement
    X = np.zeros([dim,dim])
    for i in range(dim):
        for j in range(dim):
            X[i,j] = float(i+j+1)   


    # get matrix Y with all the results 
    Y = np.zeros([dim,dim])
    for i, a in enumerate(dCarpet[keys0[0]][keys1[0]]):
        dCase[keys0[0]][keys1[0]] = a
        for j, b in enumerate(dCarpet[keys0[1]][keys1[1]]):
            dCase[keys0[1]][keys1[1]] = b
            Y[i,j] = fun(dCase)
            #print(i, j, a,b, Y[i,j])
    print('Carpet calculation finished! Returning matrices...')
                
    return X, Y

def calcSensitivitiesOAT(fun, dCase, dSample):
    # Framework for One-at-a-Time sampling sensitivity analysis
    dResult = {}
    for k0, v0 in dSample.items():
        dResult[k0] = {}
        #print('# Variation of', k0, 'parameters')
        for k1, v1, in v0.items():
            #print('## Variation of', k1)
            listRes = []
            for i, value in enumerate(v1):
                baseline = dCase[k0][k1]
                dCase[k0][k1] = value
                listRes.append(fun(dCase))
                dCase[k0][k1] = baseline
                #if dCase['fuel']['fueltype'] == 'battery': print(listRes[-1])
            dResult[k0][k1] = listRes                        
    return dResult