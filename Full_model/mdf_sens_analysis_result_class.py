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
        
        #save if a compound was varied or not and what the values were
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
        
        if all(self._S_netR[0] == self._S_netR[1]):
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
        
        all_entries_equal = all(int(element*1e4) == int(self._rNADH[0]*1e4) for element in self._rNADH)
        if all_entries_equal:
            self._rNADH = self._rNADH[0]
        
        if all(self._rNADPH):
            all_entries_equal = all(int(element*1e4) == int(self._rNADPH[0]*1e4) for element in self._rNADPH)
            if all_entries_equal:
                self._rNADPH = self._rNADPH[0]
        else:
            self._rNADPH = self._rNADPH[0]
        
        if all(self._rFd):
            all_entries_equal = all(int(element*1e3) == int(self._rFd[0]*1e3) for element in self._rFd)
            if all_entries_equal:
                self._rFd = self._rFd[0]
        else:
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
                if not (self._compounds[j] == 'Pi' or self._compounds[j] == 'H2O'):
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
        
        conditions, var = self.get_conditions()
        
        
        #concentrations
        comp_exceptions = ['H+', 'H2O', 'rFd', 'rNADH', 'rNADPH', 'rATP', 'H2', 'CO2']
        remove = []
        for i, comp in enumerate(self._compounds):
            if comp in comp_exceptions:
                remove += [i]

        plot_compounds = np.delete(np.array(self._compounds), remove)
        
        #make figure
        fig = plt.figure(figsize=(15,10))
        
        #inidividual deltaG values of reactions
        ax1 = fig.add_subplot(211)
        #ax1.set_title('Reaction energies', fontweight='bold', va = 'bottom')
        
        #compound concentrations
        ax2 = fig.add_subplot(212)
        #ax2.set_title('Compound concentrations', fontweight='bold')
        
        for i in range(len(self._opt_conc)):
            
            
            lab1 = '\n $\Delta$G_r = '+ f'{self._totaldG[i]:.2f} kJ/mol \n MDF = {max(self._dg_prime_opt[i]):.2f} kJ/mol \n'
            if var[0] == 'Stoichiometry':
                lab2 = f'\n{var[0]} {var[1][i]} \n'
            else:
                lab2 = f'\n{var[0]} = {var[1][i]:.2e} {var[2]} \n'

            
            ax1.plot(self._reactions, self._dg_prime_opt[i], 'o', c= colours[i], label = lab1)
            
            plot_conc = np.delete(np.array(self._opt_conc[i]), remove)
            ax2.plot(plot_compounds, plot_conc, 'o',c = colours[i], label = lab2)
            ax2.plot(plot_compounds, plot_conc, '--', c = colours[i], alpha =0.4)
        
        if len(self._netreaction_eq) == 1:
            neteq = self._netreaction_eq[0] 
            plt.suptitle(neteq + '\n\n\n' + conditions, ha='center' )
        else:
            plt.suptitle(conditions, ha='center' )
            j = np.arange(0,len(self._netreaction_eq))[::-1]
            for i in range(len(self._netreaction_eq)):
                ax1.text(0.5, 1.05+(j[i]/100)+(j[i]/10), f'Stoichiometry {i+1}: {self._netreaction_eq[i]}', size=12, ha="center", c = colours[i], transform=ax1.transAxes)
        
        if self._change_resulting_from_var:
            j = np.arange(0,len(self._change_resulting_from_var[1]))[::-1]
            for i in range(len(self._change_resulting_from_var[1])):
                ax1.text(0.5, 1.05+(j[i]/100)+(j[i]/10), f'{self._change_resulting_from_var[0]} = {self._change_resulting_from_var[1][i]:.2e}', size=12, ha="center", c = colours[i], transform=ax1.transAxes)
        
            
            
        ax1.grid(True)
        ax1.xaxis.set_ticks(self._reactions)
        ax1.set_xticklabels(self._reactions, fontsize=8, rotation=90)
        ax1.set_ylabel('$\Delta$G [kJ/mol]')
        ax1.set_xlabel('Reactions')
        ax1.legend(bbox_to_anchor = (1.05, 1.04))
        
        ax2.grid()
        ax2.set_ylabel('Concentration [M]')
        ax2.set_xlabel('Compounds')
        ax2.xaxis.set_ticks(plot_compounds)
        ax2.set_xticklabels(plot_compounds, fontsize=8, rotation=90)
        ax2.plot([0, len(plot_conc)], [10**-6, 10**-6], 'k-')
        ax2.plot([0, len(plot_conc)], [10**-2, 10**-2], 'k-')
        ax2.set_yscale('log')
        
        ax2.legend(bbox_to_anchor = (1.05, 1.04))
        
        plt.tight_layout()
        
        return
    
    def get_conditions(self):
        #create string of reaction conditions
        conditions = ''
        
        check_len = [self._T, self._pH, self._maxCoA, self._maxPi, self._rNADH, self._rNADPH, self._rFd, self._dGatp, self._dGprime_hyd, self._dGprime_hydABC]
        options = ['T', 'pH',  'CoA-pool', 'Pi-pool',  'rNADH', 'rNADPH', 'rFd', 'dGatp','dGprime_hyd', 'dGprime_hydABC']
        units = ['K', '', 'M', 'M', '', '', '', 'kJ/mol', 'kJ/mol', 'kJ/mol']
        
        #get name of parameter that was varied
        var_i = [i for i,c in enumerate(check_len) if type(check_len[i]) == list]
        self._change_resulting_from_var = None
        
        if self._vary_comp: 
            var = [self._vary_comp, self._var_conc, 'M' ]
            
            
            for i in range(len(check_len)):
                #unless you're at the index of the condition that is being varied, add to conditions string
                if check_len[i]:
                    if type(check_len[i]) == float or type(check_len[i]) == int:
                        conditions += f'{options[i]} = {check_len[i]:.3f} {units[i]} \t'#.expandtabs(30)
                if (i+1)%5 == 0:
                    conditions += '\n\n'
                    
            if var_i:
                var_i = var_i[0]
                self._change_resulting_from_var = [ options[var_i], check_len[var_i], units[var_i]]
                    
        elif len(var_i) != 0:
            var_i = var_i[0]
            var = [ options[var_i], check_len[var_i], units[var_i]]
        
            #loop through all conditions
            for i in range(len(check_len)):
                #unless you're at the index of the condition that is being varied, add to conditions string
                if i != var_i and check_len[i]:
                    conditions += f'\t {options[i]} = {check_len[i]:.3f} {units[i]} '.expandtabs()
                if (i+1)%5 == 0:
                    conditions += '\n\n'
        else:
            for i in range(len(check_len)):
                #unless you're at the index of the condition that is being varied, add to conditions string
                if check_len[i]:
                    conditions += f'\t {options[i]} = {check_len[i]:.3f} {units[i]} '.expandtabs()
                if (i+1)%5 == 0:
                    conditions += '\n\n'
                    
            if len(self._S_netR) > 1:
                stoichs = [ i+1 for i in range(len(self._S_netR))]
                var = ['Stoichiometry', stoichs , '']
            

            
                
        
        comp_exceptions = ['H2O', 'rNADH', 'rNADPH', 'rFd', 'CO2', 'H2']
        for i, c in enumerate(comp_exceptions):
            if c not in options and c != 'H2O' and c != var[0]:
                i_c = self._compounds.index(c)
                conditions += f'\t {c} = {self._opt_conc[0][i_c]:.2e} M '.expandtabs()
                
        return conditions, var