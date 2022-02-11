# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 15:13:14 2021

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
glu_to_product = MDF_Analysis('glu_to_acetic_butyric_all_options') 

#%% First check the standard scenario with and without physiological boundaries
glu_to_product.to_default()
result = glu_to_product.execute_mdf_basis(set_fixed_c=True, phys_bounds = True)
result.plot_results()
conc_val, dg_val = result.results_table()

#%% no phys. bounds
# result = glu_to_product.execute_mdf_basis(set_fixed_c=True, phys_bounds = False)
# result.plot_results()

#%%
var_Pi = [30e-3, 40e-3, 50e-3]

result = glu_to_product.sensitivity_analysis('Pi-pool', var_Pi)
conditions, var = result.get_conditions()
result.plot_results()

#%%
var_dG_hydABC = [0, -2, -4, -6, -10]
#glu_to_product.set_maxPi(50e-3)
result = glu_to_product.sensitivity_analysis('dGprime_hydABC', var_dG_hydABC)
conditions, var = result.get_conditions()
result.plot_results()


#%%
var_rFd = [0.01, 0.1, 0.5, 1, 2, 5, 10]

result = glu_to_product.sensitivity_analysis('rFd', var_rFd)
conditions, var = result.get_conditions()
result.plot_results()

#%%
# var_T = [298.15, 298.15+15, 298.15+30, 298.15+45]

# result = glu_to_product.sensitivity_analysis('T', var_T)
# conditions, var = result.get_conditions()
# result.plot_results()

#%%
var_cH2 = [1e-6, 1e-4, 1e-2]

result = glu_to_product.sensitivity_analysis('H2', var_cH2, vary_compound_conc = True)
conditions, var = result.get_conditions()
result.plot_results()



#%% Try different stoichiometries, run per cell after updating excel sheet
glu_to_product = MDF_Analysis('glu_to_acetic_butyric_all_options') 
res1 = glu_to_product.execute_mdf_basis(set_fixed_c=True, phys_bounds = True)

#%%
glu_to_product = MDF_Analysis('glu_to_acetic_butyric_all_options') 
res2 = glu_to_product.execute_mdf_basis(set_fixed_c=True, phys_bounds = True)

#%%
glu_to_product = MDF_Analysis('glu_to_acetic_butyric_all_options') 
res3 = glu_to_product.execute_mdf_basis(set_fixed_c=True, phys_bounds = True)

#%%
res = [res1, res2, res3]

check = MDF_Sens_Analysis_Result(res)
check.plot_results()
