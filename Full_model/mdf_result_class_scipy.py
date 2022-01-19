# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 16:15:38 2021

@author: marit
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from datafile import R

#%%

class MDF_Result(object):
    
    def __init__(self, opt_conc, dg_prime_opt, overall_dg_prime, dG0_path, reactions, compounds, fixed_c_names, S_netR, rATP_in_reaction, 
                 T, pH, ph2, pCO2, maxCoA, maxPi, rNADH, rNADPH, rFd, dGatp, dGatp0, ATP_yield):
        
        self._opt_conc = opt_conc
        self._dg_prime_opt = dg_prime_opt
        self._totaldG = overall_dg_prime
        self._dg0 = dG0_path
        self._reactions = reactions
        self._compounds = compounds
        self._fixed_c_names = fixed_c_names
        self._S_netR = S_netR
        self._rATP_in_reaction = rATP_in_reaction
        self._T = T
        self._pH = pH
        self._ph2 = ph2
        self._p_co2 = pCO2
        self._maxCoA = maxCoA
        self._maxPi = maxPi
        self._rNADH = rNADH
        self._rNADPH = rNADPH
        self._rFd = rFd
        self._dGatp = dGatp
        self._dGatp0 = dGatp0
        self._yATP = ATP_yield
        
        
        
        i_Pi = self._compounds.index('Pi')
        cH = 10**-self._pH
        self._rATP = np.exp( (self._dGatp - self._dGatp0) / (R*self._T)) * self._opt_conc[i_Pi] * (cH)
        
        #create string of net reaction to add to figure title
        sub = ''
        prod = ''
        for j in range(len(self._S_netR)):
            if not self._compounds[j] == 'Pi':
                 #to account for round-off errors: values to are ~0 should not be included in the reaction
                 if abs(self._S_netR[j]) > 1e-9:
                     if self._S_netR[j] < 0:
                         sub += f'{abs(self._S_netR[j]):.2f} {self._compounds[j]} + '
                     if self._S_netR[j] > 0:
                        prod += f' {abs(self._S_netR[j]):.2f} {self._compounds[j]} + '
            
        
        if self._yATP > 0:
            prod += f' {abs(self._yATP):.2f} ATP +'
        else:
            sub += f' {abs(self._yATP):.2f} ATP +'
        
        sub = sub[0:-2]   #remove extra plus at the end for both sides of the reaction
        prod = prod[0:-2]
        
        #save the netto reaction equation as attribute of object
        self._netreaction_eq = sub + u'\u279E' + prod
        
        self._key_info = [self._netreaction_eq, compounds, reactions, dG0_path]
        
        return

    def plot_results(self):
        """     Function that plots results in 1 figure.
                Shows netto reaction equation and reaction conditions.
                2 subplots containing the deltaG values of the reactions and the optimized concentrations."""
        
        
        #create string of reaction conditions
        conditions = f'pH = {self._pH}'.expandtabs()
        comp_exceptions = ['H2O', 'rNADH', 'rNADPH', 'rFd', 'CO2', 'H2']
        
        for comp in comp_exceptions:
            if comp in self._compounds and comp != 'H2O':
                if comp[0] == 'r':
                    i = self._compounds.index(comp)
                    conditions += f'\t {comp} = {self._opt_conc[i]:.3f}'.expandtabs()
                else:
                    i = self._compounds.index(comp)
                    conditions += f'\t [{comp}] = {self._opt_conc[i]:.2e} M'.expandtabs()

        conditions += f'\t Pi-pool = {self._maxPi:.2f} M \t CoA-pool = {self._maxCoA:.2f} M'.expandtabs()
        
        #make figure
        fig = plt.figure(figsize=(13,9))
        plt.suptitle(self._netreaction_eq + '\n\n' + 
                     f'Overall reaction energy: {self._totaldG:.2f} kJ/mol' + 
                     f'\n MDF = {max(self._dg_prime_opt):.2f} kJ/mol'  +
                     '\n\n' + conditions)

        #inidividual deltaG values of reactions
        ax2 = fig.add_subplot(211)
        ax2.title.set_text('Reaction energies')
        ax2.plot(self._reactions, self._dg0, 'rs', label='dG0 [kJ/mol]', alpha=0.4)
        ax2.plot(self._reactions, self._dg_prime_opt, 'bo', label = "dG' optimized [kJ/mol]")
        ax2.grid()
        ax2.set_ylabel('$\Delta$G [kJ/mol]')
        ax2.set_xlabel('Reactions')
        ax2.xaxis.set_ticks(self._reactions)
        ax2.set_xticklabels(self._reactions, fontsize=8, rotation=90)
        ax2.legend()

        #concentrations
        remove = []
        for i, comp in enumerate(self._compounds):
            if comp in comp_exceptions:
                remove += [i]

        plot_compounds = np.delete(np.array(self._compounds), remove)
        plot_conc = np.delete(np.array(self._opt_conc), remove)
                
        ax3 = fig.add_subplot(212)
        ax3.title.set_text('Compound concentrations')
        
        col = ['b']*len(plot_compounds)
        
        for i, c in enumerate(plot_compounds):
            if c in self._fixed_c_names:
                col[i] = 'r'
        
        for i in range(len(plot_compounds)):
            ax3.scatter(plot_compounds[i], plot_conc[i], c=col[i])
            
        ax3.plot(plot_compounds, plot_conc, 'b--', alpha =0.4)
        ax3.grid()
        ax3.set_ylabel('Concentration [M]')
        ax3.set_xlabel('Compounds')
        ax3.xaxis.set_ticks(plot_compounds)
        ax3.set_xticklabels(plot_compounds, fontsize=8, rotation=90)
        ax3.plot([0, len(plot_conc)], [10**-6, 10**-6], 'k-')
        ax3.plot([0, len(plot_conc)], [10**-2, 10**-2], 'k-')
        ax3.set_yscale('log')
        
        legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor='b', label='Opt. conc.', markersize=7),
                           Line2D([0], [0], marker='o', color='w', markerfacecolor='r', label='Fixed conc.', markersize=7)]

        ax3.legend(handles=legend_elements, loc='upper right')

        plt.tight_layout()
        
        return fig
    
    def results_table(self):
        dict_conc = {}
        for i, comp in enumerate(self._compounds):
            dict_conc[comp] = f'{self._opt_conc[i]:.3e}'
            
        dict_dg = {}
        for i, reac in enumerate(self._reactions):
            dict_dg[reac] = f'{self._dg_prime_opt[i]:.3f}'
        
        return dict_conc, dict_dg
    
    def fluxforce_efficiacy(self):
        """     Calculates flux force efficiacy of a reaction.
                The closer the value to 1, the less back-flux there is."""
        FFE = ( np.exp(-self._dg_prime_opt / (R * self._T) ) - 1 ) / ( np.exp(-self._dg_prime_opt / (R * self._T) ) + 1 )
        
        dict_FFE = {}
        for i, reac in enumerate(self._reactions):
            dict_FFE[reac] = f'{FFE[i]:.3f}'
        
        return dict_FFE
    
    def export_results(self):
        
        return