# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 12:51:00 2021

@author: marit
"""
#%%
import sys
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import root

#if MEPfunctions.py and datafile.py not in same folder, add the required folder to path
sys.path.append('C:\\Users\marit\Documents\LST\MSc\MEP\Scipy MDF\MDF-ECM')
sys.path.append('C:/Users/marit/Documents/LST/MSc/MEP/Scipy MDF/MDF-ECM/Full_model')

#%%
from datafile import (
    def_c_max,
    def_c_min,
    F,
    R)

from pathway_class_scipy_cc import Pathway_cc
import warnings

#%%
class Equilibrium_Analysis(Pathway_cc):
    

    
    def eq_model(self):
        
        def find_conc(ln_conc):
            i_Pi = self._compounds.index('Pi')
            cH = 10**-self._pH
            
            rATP = np.exp( (self._dGatp - self._dGatp0) / (R*self._T)) * np.exp(ln_conc[i_Pi]) * (cH)
            
            #set dg_prime as a function of dg0_prime and variables ln_conc
            dg_prime = self._dg0 + ( R*self._T * self._stoich.T @ ln_conc ) + ( R * self._T * self._rATP_in_reaction * np.log(rATP) )
             
            #totalPi = total amount of compounds carrying phosphate
            #maxPi = phosphate pool
            
            #get indices of compounds carrying phosphate
            i_p = np.where(self._element_comp[:,-2] != 0)[0]
                    
            totalPi = 0
            for i in i_p:
                #don't include ATP in phosphate constraint
                if not self._compounds[i] == 'rATP]':
                    #x_p is the amount of phosphate a compound carries
                    x_p = self._element_comp[i,-2]
                    #multiply concentration by amount of phosphate in molecule
                    compconc = x_p * np.exp(ln_conc[i])
                    #add moles of phosphate together
                    totalPi += compconc
                    
            #totalCoA= total amount of compounds carrying CoA
            #maxCoA = CoA pool
            
            #get indices of compounds carrying CoA
            i_coA = []
            for i, string in enumerate(self._compounds): 
                if 'CoA' in string: 
                    i_coA += [i]
                    
            totalCoA = 0
            #add all concentrations together
            for i in i_coA:
                compconc = np.exp(ln_conc[i])
                totalCoA += compconc
           
            #constraints that should equal to zero: total pool = total concentration of compounds carrying the element
            consPi = totalPi - self._maxPi
            consCoA = totalCoA - self._maxCoA
            
            #create single array with equations that should all be equal to 0
            y = [dg_prime, consPi, consCoA]
            
            return y
        
        conc0 = [np.log(1e-4)] * len(self._compounds)
        solve = root(find_conc, conc0, method='lm')

        self._eq_result_conc = np.exp(solve.x)

        self._eq_result_dg = self._dg0 + ( R*self._T * self._stoich.T @ solve.x )
        
        return self._eq_result_conc, self._eq_result_dg
    
glu_to_ac_ebif = Equilibrium_Analysis('try') 

# dg0 = glu_to_ac_ebif._dg0
# T = glu_to_ac_ebif._T
# stoich = glu_to_ac_ebif._stoich
# rATP_in_reaction = glu_to_ac_ebif._rATP_in_reaction
# element_comp = glu_to_ac_ebif._element_comp
# compounds = glu_to_ac_ebif._compounds

# maxPi = glu_to_ac_ebif._maxPi
# maxCoA = glu_to_ac_ebif._maxCoA

# # ln_conc = np.zeros(len(glu_to_ac_ebif._compounds))
# # y = dg0 + ( R*T * stoich.T @ ln_conc )

# def find_conc(ln_conc):
    
#     y = dg0 + ( R*T * stoich.T @ ln_conc ) #+ ( R * T * rATP_in_reaction * np.log(rATP) )
    
    
#     #totalPi = total amount of compounds carrying phosphate
#     #maxPi = phosphate pool
    
#     #get indices of compounds carrying phosphate
#     i_p = np.where(element_comp[:,-2] != 0)[0]
            
#     totalPi = 0
#     for i in i_p:
#         #don't include ATP in phosphate constraint
#         if not compounds[i] == 'rATP]':
#             #x_p is the amount of phosphate a compound carries
#             x_p = element_comp[i,-2]
#             #multiply concentration by amount of phosphate in molecule
#             compconc = x_p * np.exp(ln_conc[i])
#             #add moles of phosphate together
#             totalPi += compconc
    
#     #totalCoA= total amount of compounds carrying CoA
#     #maxCoA = CoA pool
    
#     #get indices of compounds carrying CoA
#     i_coA = []
#     for i, string in enumerate(compounds): 
#         if 'CoA' in string: 
#             i_coA += [i]
            
#     totalCoA = 0
#     #add all concentrations together
#     for i in i_coA:
#         compconc = np.exp(ln_conc[i])
#         totalCoA += compconc
    
#     consPi = totalPi - maxPi
#     consCoA = totalCoA - maxCoA
    
#     y = np.append(y, [consPi, consCoA])
    
#     return y

# #%%
# conc0 = [np.log(1e-5)] * len(compounds)
# solve = root(find_conc, conc0, method='lm')

# res = np.exp(solve.x)

# y = dg0 + ( R*T * stoich.T @ solve.x )

#%% 
#glu_to_ac_ebif = Equilibrium_Analysis('try') 
test = glu_to_ac_ebif.eq_model()