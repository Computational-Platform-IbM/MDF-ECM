# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 18:28:25 2021

@author: marit
"""
#%%
from mdf_class_scipy import MDF_Analysis
import matplotlib.pyplot as plt
import numpy as np

#%%

glu_to_ac = MDF_Analysis('glu_to_acetate_ratios')
glu_to_ac_ebif = MDF_Analysis('glu_to_ac_e-bfc_ratios')

glu_to_ac.set_solver_tol(1e-7)
glu_to_ac_ebif.set_solver_tol(1e-7)

result = glu_to_ac.mdf_fixed_conc(fixed_rNADH=False)

rNADH = result._rNADH
glu_to_ac_ebif.set_rNADH(rNADH)

result_ebif = glu_to_ac_ebif.mdf_fixed_conc()

result.plot_results()
result_ebif.plot_results()

varpH2 = [1e-2, 1, 5]

for pH2 in varpH2:
    glu_to_ac_ebif.set_p_h2(pH2)
    
    result = glu_to_ac_ebif.mdf_fixed_conc()
    #result.plot_results()
    
glu_to_ac_ebif.set_rNADH(0.05)
for pH2 in varpH2:
    glu_to_ac_ebif.set_p_h2(pH2)
    
    result = glu_to_ac_ebif.mdf_fixed_conc()
    #result.plot_results()

#%% Vary pi pool and NADH ratio
varPi = [10e-3, 20e-3, 30e-3]
var_rNADH = [5e-4, 5e-2, 0.5]

glu_to_ac.set_solver_tol(1e-7)
glu_to_ac_ebif.set_solver_tol(1e-7)

glu_to_ac.set_p_h2(1), glu_to_ac_ebif.set_p_h2(1)

conc = []
conc_ebif = []

dg_opt = []
dg_opt_ebif = []

for valPi in varPi:
    glu_to_ac.set_maxPi(valPi)
    glu_to_ac_ebif.set_maxPi(valPi)
    for rNADH in var_rNADH:
        glu_to_ac.set_rNADH(rNADH)
        glu_to_ac_ebif.set_rNADH(rNADH)
        
        result = glu_to_ac.mdf_fixed_conc()
        conc_val, dg_val = result.results_table()
        conc += [{'maxPi' : valPi, 'rNADH': rNADH, 'conc': conc_val}]
        dg_opt += [{'maxPi' : valPi, 'rNADH': rNADH, 'dg': dg_val}]
        
        result_ebif = glu_to_ac_ebif.mdf_fixed_conc()
        conc_val, dg_val = result_ebif.results_table()
        conc_ebif += [{'maxPi' : valPi, 'rNADH': rNADH, 'conc': conc_val}]
        dg_opt_ebif += [{'maxPi' : valPi, 'rNADH': rNADH, 'dg': dg_val}]

#%%
#make figure
plt.figure(figsize=(13,9))
plt.suptitle(result._netreaction_eq)

#inidividual deltaG values of reactions
plt.title('Optimized reaction energies')

for i in range(0,len(conc)):
    plot_conc = list(conc[i]['conc'].values())
    plot_dg = []
    for val in list(dg_opt[i]['dg'].values()):
        plot_dg += [float(val)]
    
    plt.plot(result._reactions, plot_dg, 'o', label = f"rNADH = {conc[i]['rNADH']}, maxPi = {conc[i]['maxPi']}")

plt.grid()
plt.ylabel('$\Delta$G [kJ/mol]')
plt.xlabel('Reactions')
plt.ylim([-20, 10])
plt.legend()
