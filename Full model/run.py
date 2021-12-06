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

glu_to_ac = MDF_Analysis('try') #'glu_to_ac_e-bfc_ratios_new') #'glu_to_acetate_ratios_new') #

#%%

result = glu_to_ac.execute_mdf_basis(set_fixed_c=True, fixed_rNADH = True, phys_bounds = True)
result.plot_results()
print(result._dg_prime_opt)
print(sum(result._dg_prime_opt))
print(result._rATP)

result = glu_to_ac.execute_mdf_basis(set_fixed_c=True, fixed_rNADH = True, phys_bounds = False)
result.plot_results()
print(result._dg_prime_opt)
print(sum(result._dg_prime_opt))
print(result._rATP)

result = glu_to_ac.execute_mdf_basis(set_fixed_c=True, fixed_rNADH = False, phys_bounds = True)
result.plot_results()
print(result._dg_prime_opt)
print(sum(result._dg_prime_opt))
print(result._rATP)


#%% Vary NADH ratio
#var_rNADH = [3e-6, 1/240, 5e-2, 0.16]
var_rNADH = [0.4, 0.05]
var_Pi = [10e-3, 20e-3, 30e-3, 50e-3]

conc = []

dg_opt = []



for val in var_Pi:
    glu_to_ac.set_maxPi(val)
    
    for val2 in var_rNADH:
        glu_to_ac.set_rNADH(val2)
    
        result = glu_to_ac.execute_mdf_basis(set_fixed_c=True, fixed_rNADH = True, phys_bounds = True)
        conc_val, dg_val = result.results_table()
        conc += [{'maxPi': val, 'rNADH' : val2, 'conc': conc_val}]
        dg_opt += [{'maxPi': val, 'rNADH' : val2, 'dg': dg_val}]

#%%
#make figure
fig = plt.figure(figsize=(13,9))
plt.suptitle(result._netreaction_eq)

#inidividual deltaG values of reactions
ax2 = fig.add_subplot(211)
ax3 = fig.add_subplot(212)

#inidividual deltaG values of reactions
plt.title('Optimized reaction energies')

#concentrations
comp_exceptions = ['H+', 'H2O', 'H2', 'rFd', 'rNADH', 'rNADPH', 'rATP']
remove = []
for i, comp in enumerate(result._compounds):
    if comp in comp_exceptions:
        remove += [i]

plot_compounds = np.delete(np.array(result._compounds), remove)

for i in range(0,len(conc)):
    plot_conc = list(conc[i]['conc'].values())
    plot_conc = [float(val) for val in plot_conc]
    plot_conc = np.delete(np.array(plot_conc), remove)
    
    plot_dg = [float(val) for val in list(dg_opt[i]['dg'].values())]
    
    ax2.plot(result._reactions, plot_dg, 'o', label = f"maxPi = {conc[i]['maxPi']:.2e}, rNADH = {conc[i]['rNADH']:.2e}")
    ax3.plot(plot_compounds, plot_conc, 'o--', label = f"maxPi = {conc[i]['maxPi']:.2e}, rNADH = {conc[i]['rNADH']:.2e}")

ax2.grid(True)
ax3.grid(True)
ax2.set_ylabel('$\Delta$G [kJ/mol]')
ax2.set_xlabel('Reactions')
#ax2.set_ylim([-40, 5])
ax2.legend()

ax3.set_xticklabels(plot_compounds, fontsize=8, rotation=90)
ax3.plot([0, len(plot_conc)], [10**-6, 10**-6], 'k-')
ax3.plot([0, len(plot_conc)], [10**-2, 10**-2], 'k-')
ax3.set_yscale('log')
ax3.legend()

plt.tight_layout()
