# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 13:48:41 2021

@author: marit
"""

#%%
import numpy as np
import pandas as pd
import os
#import imageio
from typing import Dict, List
import matplotlib.pyplot as plt
#import tqdm

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

    rel_flux = data_ar[-1,8:]
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

#%% Functions to save plots and create gifs - adapted from Chiel

def animate_optimization(i: int, xlim: List[float], ylim: List[float], bac: Dict):
    """Create and save a figure of the bacteria at a certain point in time.
    Args:
        i (int): [description]
        xlim (List[float]): [description]
        ylim (List[float]): [description]
    Returns:
        string: name of the file in which the figure is saved.
    """

    # filename = f'{directory}/{i}.png'
    # plt.savefig(filename)
    # plt.close()
    
    #generate_gif(args)

    return #filename


# %%
def generate_gif(args: Dict):
    
    return
    
#%% 

def comparison_plot(dgs: Dict, conc: Dict, mdf: List, total_dG: List, key_info: List):
    netreaction_eq, compounds, reactions, dg0 = key_info
    
    colours = ['mediumblue', 'darkorchid', 'palevioletred', 
               'firebrick', 'tan', 'darkcyan', 'seagreen']
    
    #make figure
    fig = plt.figure(figsize=(13,9))
    plt.suptitle(netreaction_eq)

    #inidividual deltaG values of reactions
    ax2 = fig.add_subplot(211)
    ax2.title.set_text('Reaction energies')
    #ax2.plot(reactions, dg0, 'o', label = "$\Delta$G_0")
    ax3 = fig.add_subplot(212)


    #concentrations
    comp_exceptions = ['H+', 'H2O', 'rFd', 'rNADH', 'rNADPH', 'rATP']
    remove = []
    for i, comp in enumerate(compounds):
        if comp in comp_exceptions:
            remove += [i]

    plot_compounds = np.delete(np.array(compounds), remove)

    for i in range(0,len(conc)):
        plot_conc = list(conc[i]['conc'].values())
        plot_conc = [float(val) for val in plot_conc]
        plot_conc = np.delete(np.array(plot_conc), remove)
        
        plot_dg = [float(val) for val in list(dgs[i]['dg'].values())]
        
        lab1 = f'Total $\Delta$G = {total_dG[i]:.2f} kJ/mol'
        lab2 = ''
        
        for key in conc[i]:
            if key != 'conc':
                lab2 += f'\n {key} = {conc[i][key]:.2e}, '
        
        lab2 += f'\n MDF = {mdf[i]:.2f} kJ/mol'
        
        ax2.plot(reactions, plot_dg, 'o', color = colours[i], label = lab1)
        ax3.plot(plot_compounds, plot_conc, 'o', color = colours[i], label = lab2)
        ax3.plot(plot_compounds, plot_conc, '--', color = colours[i],  alpha = 0.4)

    ax2.grid(True)
    ax3.grid(True)
    ax2.set_ylabel('$\Delta$G [kJ/mol]')
    ax2.set_xlabel('Reactions')
    ax2.legend()
    
    ax3.title.set_text('Compound concentrations')
    ax3.xaxis.set_ticks(plot_compounds)
    ax3.set_xticklabels(plot_compounds, fontsize=8, rotation=90)
    ax3.plot([0, len(plot_conc)], [10**-6, 10**-6], 'k-')
    ax3.plot([0, len(plot_conc)], [10**-2, 10**-2], 'k-')
    ax3.set_yscale('log')
    ax3.legend()

    plt.tight_layout()
    
    return
