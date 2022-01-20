# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 10:33:38 2022

@author: marit
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from datafile import R

#%%

class MDF_Sens_Analysis_Result(object):
    
    def __init__(self, result_objects, vary_comp_conc = None, var_conc = None):
        
        self._opt_conc          = []
        self._dg_prime_opt      = []
        self._totaldG           = []
        self._reactions         = []
        self._compounds         = []
        #self._fixed_c_names     = []
        self._S_netR            = []
        self._rATP_in_reaction  = []
        self._T                 = []
        self._pH                = []
        #self._ph2               = []
        #self._p_co2             = []
        self._maxCoA            = []
        self._maxPi             = []
        self._rNADH             = []
        self._rNADPH            = []
        self._rFd               = []
        self._dGatp             = []
        self._dGatp0            = []
        self._rATP              = []
        self._yATP              = []
        self._dGprime_hyd       = []
        self._dGprime_hydABC    = []
        
        if vary_comp_conc:
            self._vary_comp = vary_comp_conc
            self._var_conc = var_conc
        
        
        for i in range(len(result_objects)):
            self._opt_conc          += [result_objects[i]._opt_conc]
            self._dg_prime_opt      += [result_objects[i]._dg_prime_opt]
            self._totaldG           += [result_objects[i]._totaldG]
            self._reactions         += [result_objects[i]._reactions]
            self._compounds         += [result_objects[i]._compounds]
            #self._fixed_c_names     += [result_objects[i]._fixed_c_names]
            self._S_netR            += [result_objects[i]._S_netR]
            self._yATP              += [result_objects[i]._yATP ]
            self._rATP_in_reaction  += [result_objects[i]._rATP_in_reaction]
            
            #variables
            self._T                 += [result_objects[i]._T]
            self._pH                += [result_objects[i]._pH]
            #self._ph2               += [result_objects[i]._ph2]
            #self._p_co2             += [result_objects[i]._p_co2 ]
            self._maxCoA            += [result_objects[i]._maxCoA ]
            self._maxPi             += [result_objects[i]._maxPi ]
            self._rNADH             += [result_objects[i]._rNADH]
            self._rNADPH            += [result_objects[i]._rNADPH]
            self._rFd               += [result_objects[i]._rFd]
            self._dGatp             += [result_objects[i]._dGatp]
            self._dGatp0            += [result_objects[i]._dGatp0]
            
            i_Pi = result_objects[i]._compounds.index('Pi')
            cH = 10**-result_objects[i]._pH
            self._rATP              += [np.exp( (result_objects[i]._dGatp - result_objects[i]._dGatp0) / (R*result_objects[i]._T)) * result_objects[i]._opt_conc[i_Pi] * (cH)]
            
            self._dGprime_hyd       += [result_objects[i]._dGprime_hyd]
            self._dGprime_hydABC    += [result_objects[i]._dGprime_hydABC]
            
        
        #can only do plots of sensitivity analysis of pathways that have the same reactions
        #not necessarily the same stoichiometry!
        
        if not self._reactions[0] == self._reactions[1]:
            raise ValueError('Pathways of the analyses are not compatible!')
        else:
            self._reactions = self._reactions[0]
            self._compounds = self._compounds[0]
        
        if all(self._S_netR[0]) == all(self._S_netR[1]):
            self._S_netR = [self._S_netR[0]]
            
        all_entries_equal = all(element == self._T[0] for element in self._T)
        if all_entries_equal:
            self._T = self._T[0]
        
        all_entries_equal = all(element == self._pH[0] for element in self._pH)
        if all_entries_equal:
            self._pH = self._pH[0]
        
        all_entries_equal = all(element == self._maxPi[0] for element in self._maxPi)
        if all_entries_equal:
            self._maxPi = self._maxPi[0]
        
        all_entries_equal = all(element == self._maxCoA[0] for element in self._maxCoA)
        if all_entries_equal:
            self._maxCoA = self._maxCoA[0]
        
        all_entries_equal = all(int(element*1e6) == int(self._rNADH[0]*1e6) for element in self._rNADH)
        if all_entries_equal:
            self._rNADH = self._rNADH[0]
        
        all_entries_equal = all(int(element*1e6) == int(self._rNADPH[0]*1e6) for element in self._rNADPH)
        if all_entries_equal:
            self._rNADPH = self._rNADPH[0]
        
        all_entries_equal = all(int(element*1e6) == int(self._rFd[0]*1e6) for element in self._rFd)
        if all_entries_equal:
            self._rFd = self._rFd[0]
        
        all_entries_equal = all(element == self._dGatp[0] for element in self._dGatp)
        if all_entries_equal:
            self._dGatp = self._dGatp[0]
            
        all_entries_equal = all(element == self._dGatp0[0] for element in self._dGatp0)
        if all_entries_equal:
            self._dGatp0 = self._dGatp0[0]
        
        all_entries_equal = all(element == self._dGprime_hyd[0] for element in self._dGprime_hyd)
        if all_entries_equal:
            self._dGprime_hyd = self._dGprime_hyd[0]
        
        all_entries_equal = all(element == self._dGprime_hydABC[0] for element in self._dGprime_hydABC)
        if all_entries_equal:
            self._dGprime_hydABC = self._dGprime_hydABC[0]
            
        self.get_reaction_eqs()
        
        # self._key_info = [self._netreaction_eq, compounds, reactions, dG0_path]
        
        return
    
    def get_reaction_eqs(self):
        
        self._netreaction_eq = []
        #create string of net reaction to add to figure title
        for i in range(len(self._S_netR)):
            
            sub = ''
            prod = ''
            for j in range(len(self._S_netR[i])):
                if not self._compounds[j] == 'Pi':
                      #to account for round-off errors: values to are ~0 should not be included in the reaction
                      if abs(self._S_netR[i][j]) > 1e-9:
                          if self._S_netR[i][j] < 0:
                              sub += f'{abs(self._S_netR[i][j]):.2f} {self._compounds[j]} + '
                          if self._S_netR[i][j] > 0:
                            prod += f' {abs(self._S_netR[i][j]):.2f} {self._compounds[j]} + '
                
            
            if self._yATP[i] > 0:
                prod += f' {abs(self._yATP[i]):.2f} ATP +'
            else:
                sub += f' {abs(self._yATP[i]):.2f} ATP +'
            
            sub = sub[0:-2]   #remove extra plus at the end for both sides of the reaction
            prod = prod[0:-2]
            
            #save the netto reaction equation as attribute of object
            self._netreaction_eq += [sub + u'\u279E' + prod]
        return
    
    def plot_results(self):
        colours = ['mediumblue', 'darkorchid', 'palevioletred', 
                   'firebrick', 'tan', 'darkcyan', 'seagreen']
        
        #create string of reaction conditions
        conditions = f'pH = {self._pH}'.expandtabs()
        comp_exceptions = ['H2O', 'rNADH', 'rNADPH', 'rFd', 'CO2', 'H2']
        
        # for comp in comp_exceptions:
        #     if comp in self._compounds and comp != 'H2O':
        #         if comp[0] == 'r':
        #             i = self._compounds.index(comp)
        #             conditions += f'\t {comp} = {self._opt_conc[i]:.3f}'.expandtabs()
        #         else:
        #             i = self._compounds.index(comp)
        #             conditions += f'\t [{comp}] = {self._opt_conc[i]:.2e} M'.expandtabs()

        # conditions += f'\t Pi-pool = {self._maxPi:.2f} M \t CoA-pool = {self._maxCoA:.2f} M'.expandtabs()
        
        #make figure
        fig = plt.figure(figsize=(13,9))
        
        neteq = self._netreaction_eq[0]
        plt.suptitle(neteq )#+ '\n\n' + 
                     # f'Overall reaction energy: {self._totaldG:.2f} kJ/mol' + 
                     # f'\n MDF = {max(self._dg_prime_opt):.2f} kJ/mol'  +
                     # '\n\n' + conditions)
        
        #concentrations
        comp_exceptions = ['H+', 'H2O', 'rFd', 'rNADH', 'rNADPH', 'rATP']
        remove = []
        for i, comp in enumerate(self._compounds):
            if comp in comp_exceptions:
                remove += [i]

        plot_compounds = np.delete(np.array(self._compounds), remove)
        
        #make figure
        fig = plt.figure(figsize=(15,9))
        plt.suptitle(self._netreaction_eq)

        #inidividual deltaG values of reactions
        ax2 = fig.add_subplot(211)
        ax2.set_title('Reaction energies', fontweight='bold')
        ax3 = fig.add_subplot(212)
        ax3.set_title('Compound concentrations', fontweight='bold')
        
        
        return