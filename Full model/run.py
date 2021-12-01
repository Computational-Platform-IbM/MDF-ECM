# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 18:28:25 2021

@author: marit
"""
#%%
import sys
#if MEPfunctions.py and datafile.py not in same folder, add the required folder to path
sys.path.append('C:\\Users\marit\Documents\LST\MSc\MEP\Scipy MDF\MDF-ECM')

from mdf_class_scipy import MDF_Analysis
import matplotlib.pyplot as plt
import numpy as np

#%%

glu_to_ac = MDF_Analysis('glu_to_acetate_ratios_new')
result = glu_to_ac.mdf_fixed_conc()
result.plot_results()
result = glu_to_ac.mdf_fixed_conc(fixed_rNADH=False)
result.plot_results()

#%%

glu_to_ac_ebif = MDF_Analysis('glu_to_ac_e-bfc_ratios_new')
result_ebif = glu_to_ac_ebif.mdf_fixed_conc()
result_ebif.plot_results()
result_ebif = glu_to_ac_ebif.mdf_fixed_conc(fixed_rNADH=False)
result_ebif.plot_results()


#%% Vary pi pool and NADH ratio
varPi = [10e-3, 20e-3, 30e-3]
var_rNADH = [5e-4, 5e-2, 0.5]

#glu_to_ac_ebif.set_solver_tol(1e-7)

conc = []
conc_ebif = []

dg_opt = []
dg_opt_ebif = []

for valPi in varPi:
    glu_to_ac_ebif.set_maxPi(valPi)
    for rNADH in var_rNADH:
        glu_to_ac_ebif.set_rNADH(rNADH)
        
        result_ebif = glu_to_ac_ebif.mdf_fixed_conc()
        conc_val, dg_val = result_ebif.results_table()
        conc_ebif += [{'maxPi' : valPi, 'rNADH': rNADH, 'conc': conc_val}]
        dg_opt_ebif += [{'maxPi' : valPi, 'rNADH': rNADH, 'dg': dg_val}]

#%%
#make figure
plt.figure(figsize=(13,9))
plt.suptitle(result_ebif._netreaction_eq)

#inidividual deltaG values of reactions
plt.title('Optimized reaction energies')

for i in range(0,len(conc_ebif)):
    plot_conc = list(conc_ebif[i]['conc'].values())
    plot_dg = []
    for val in list(dg_opt_ebif[i]['dg'].values()):
        plot_dg += [float(val)]
    
    plt.plot(result_ebif._reactions, plot_dg, 'o', label = f"rNADH = {conc_ebif[i]['rNADH']}, maxPi = {conc_ebif[i]['maxPi']}")

plt.grid()
plt.ylabel('$\Delta$G [kJ/mol]')
plt.xlabel('Reactions')
plt.ylim([-20, 10])
plt.legend()
