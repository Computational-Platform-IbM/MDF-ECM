# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 13:48:41 2021

@author: marit
"""

#%%
import numpy as np
import pandas as pd
import os

#%% import pathway info
def importpath(filename):
    # Read the data from the excel sheet
    cwd = os.getcwd()

    # add file name incl extension below
    path = cwd + filename
    # specify which sheet you want to pull data from
    data = pd.read_excel(path, 'Blad1', header=0, index_col=0)

    # create np array from the dataframe
    data_ar = np.array(data)

    # create array containing just the compound identifiers
    compound_ids = data_ar[:-1,0]
    data_ar = np.delete(data_ar, 0, axis=1)

    data_ar = np.array(data_ar, dtype=np.float64)

    #make sure no empty columns are included in the matrix
    empty = True
    i_empty = 0
    while empty == True:
        if np.isnan(data_ar[:,-1]).all():
            data_ar = np.delete(data_ar, -1, axis=1)
            i_empty += -1
        else:
            empty = False

    # create list of reaction names and compounds
    reactions = list(data.columns)
    #remove empty entries from the reaction list
    if i_empty == 0:
        reactions = reactions[9:]
    else:    
        reactions = reactions[9:i_empty]
        
    compounds = list(data.index)
    compounds.pop()             #remove relative flux from compound list

    rel_flux = data_ar[-1,9::]
    data_ar = np.delete(data_ar, -1, axis=0)

    element_bal = data.iloc[:-1, 1:7]
    element_balance = np.array(element_bal)

    #fixed_c
    fixed_c = data_ar[:,6]

    #create array for stoichiometric coefficients of netto pathway reaction
    S_netR = data_ar[:,7]
    S_netR[np.isnan(S_netR)] = 0

    # create stoichiometric matrix with just the reaction coefficiets: remove first three columns 
    stoich = data_ar[:, 8:]
    stoich[np.isnan(stoich)] = 0
    
    return reactions, compounds, element_balance, fixed_c, compound_ids, S_netR, stoich, rel_flux
