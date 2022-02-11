# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 09:06:04 2022

@author: marit
"""
#%%
import sys
#if MEPfunctions.py and datafile.py not in same folder, add the required folder to path
sys.path.append('C:\\Users\marit\Documents\LST\MSc\MEP\Scipy MDF\MDF-ECM')
sys.path.append('C:/Users/marit/Documents/LST/MSc/MEP/Scipy MDF/MDF-ECM/Full_model')

from mdf_class_scipy import MDF_Analysis
from mdf_sens_analysis_result_class import MDF_Sens_Analysis_Result
from MEPfunctions import comparison_plot
import numpy as np

#%% Initialization of pathway analysis
CE = MDF_Analysis('chain_elongation_to_but_and_cap_NADPH_v2') 

#%% First check the standard scenario with and without physiological boundaries
CE.to_default()
CE.set_rNADH(0.3)
CE.set_rNADPH(1.5)
result = CE.execute_mdf_basis(set_fixed_c=True, phys_bounds = False, user_defined_rNADH = True)#, user_defined_rNADPH = True)
result.plot_results()
conc_val, dg_val = result.results_table()

#%%
var_CoA = [20e-3, 30e-3, 40e-3, 50e-3]

CE.set_rNADH(0.3)
CE.set_rNADPH(1.5)
result = CE.sensitivity_analysis('CoA-pool', var_CoA, phys_bounds=False, user_defined_rNADH = True, user_defined_rNADPH = True)
conditions, var = result.get_conditions()
result.plot_results()

#%%
var_rNADH = [0.5, 0.3, 0.1, 0.02, 0.005]

#CE.set_rNADH(0.3)
CE.set_rNADPH(1.5)
result = CE.sensitivity_analysis('rNADH', var_rNADH, user_defined_rNADH = True, phys_bounds=False, user_defined_rNADPH = True)
result.plot_results()

#%%
var_rFd = [0.01, 0.1, 1, 10]

CE.set_rNADH(0.3)
CE.set_rNADPH(1.5)
result = CE.sensitivity_analysis('rFd', var_rFd, set_fixed_c=True, phys_bounds = False, user_defined_rFd = True, user_defined_rNADH = True, user_defined_rNADPH = True)
result.plot_results()


#%%
var_rNADPH = [0.1, 1, 10, 50, 100]

CE.to_default()
CE.set_rNADH(0.3)
CE.set_rNADPH(1.5)
result = CE.sensitivity_analysis('rNADPH', var_rNADPH, user_defined_rNADPH = True,user_defined_rNADH = True)
result.plot_results()


#%%
var_cH2 = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2]

CE.to_default()
CE.set_rNADH(0.3)
CE.set_rNADPH(1.5)
result = CE.sensitivity_analysis('H2', var_cH2, vary_compound_conc = True, phys_bounds = False, user_defined_rNADH = True)
conditions, var = result.get_conditions()
result.plot_results()

#%%
# CE1 = MDF_Analysis('chain_elongation_to_but_and_cap_NADPH_v2') 
# res1 = CE1.execute_mdf_basis(set_fixed_c=True, phys_bounds = True)

# #%%
# CE2 = MDF_Analysis('chain_elongation_to_but_and_cap_NADPH_v2') 
# res2 = CE2.execute_mdf_basis(set_fixed_c=True, phys_bounds = True)

# #%%
# CE3 = MDF_Analysis('chain_elongation_to_but_and_cap_NADPH_v2') 
# res3 = CE3.execute_mdf_basis(set_fixed_c=True, phys_bounds = True)
# #%%


# #%%
# CE4 = MDF_Analysis('chain_elongation_to_but_and_cap_NADPH_v2') 
# CE4.to_default()
# res4 = CE4.execute_mdf_basis(set_fixed_c=True, phys_bounds = True)

# res4.plot_results()
# conc_val, dg_val = res4.results_table()

#%%

# res = [res2, res3, res4]

# check = MDF_Sens_Analysis_Result(res)
# check.plot_results()
