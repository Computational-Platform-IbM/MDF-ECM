# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 10:10:52 2021

@author: marit
"""

#%%
import numpy as np

from MEPfunctions import importpath 

from datafile import (
    F,
    default_T,
    default_pH,
    default_pH2,
    default_pCO2,
    dGatp0,
    default_dGatp)

import scipy.stats
from equilibrator_api import ComponentContribution, Q_
from typing import Dict, List

#%%

class Pathway_cc(object):
    
    def __init__(self, filename, T = default_T, pH = default_pH):
        """Create a Pathway object with default settings, importing the Pathway from an Excel sheet. """
        
        #initialize pathway with default values
        self._p_h   = pH  
        self._T     = T     #K
        
        
        #set default values for CoA and Pi pool
        self._maxPi     = 20e-3     #M
        self._maxCoA    = 10e-3     #M
        
        #set default values for NADH/NAD+ and NADPH/NADP+
        self._rNADH     = 0.05
        self._rNADPH    = None
        
        #set default values for ATP production
        self._dGatp0    = dGatp0
        self._dGatp     = default_dGatp
        
        self._rFd  = None
        
        #set default tolerance for solver
        self._tol_conc  = 1e-9
        
        #if excel extension not included in filename, add it to open the file
        if not filename[-5:] == '.xlsx':
            filename += '.xlsx'
        
        #import pathway
        (self._reactions, 
         self._compounds, 
         self._element_comp,
         self._fixed_c, 
         self._compound_ids, 
         self._S_netR, 
         self._stoich, 
         self._rel_flux)     = importpath('\\' + filename)
        
        #create copy of list with all compounds and stoich matrix
        self._compounds_copy = self._compounds.copy()
        self._stoich_copy = self._stoich.copy()
        
        #default values
        self._pH2   = default_pH2   #atm
        self._pCO2   = default_pCO2   #atm
            
        
        #checks if the stoichiometry closes and automatically calculates dg0 values of all reactions
        self.check_element_balance()
        self.get_dGf_prime()
        self.calc_dG0_path()
        
        #after calculating dg0 of the pathway reactions, remove rATP from the stoichiometric matrix etc
        #save involvement of ATP/ADP in separate array
        i_rATP = self._compounds.index('rATP')
        self._rATP_in_reaction = self._stoich[i_rATP,:]
        self._netATP = self._S_netR[i_rATP]
        self._S_netR_copy = self._S_netR.copy()

        #remove rATP from matrix and arrays
        self._compounds.pop(i_rATP)
        self._stoich = np.delete(self._stoich, i_rATP, axis=0)
        self._fixed_c = np.delete(self._fixed_c, i_rATP)
        self._element_comp = np.delete(self._element_comp, i_rATP, axis=0)
        self._S_netR = np.delete(self._S_netR, i_rATP)
        
        #get number of compounds and reactions in pathway
        self._Nc, self._Nr  = self._stoich.shape
        
    @property
    def p_h(self):
        """Get the pH."""
        return self._p_h

    def set_p_h(self, value):
        """Set the pH."""
        self._p_h = value
        #if the pH is changed: update dGf' values and dg0 of reactions
        self.get_dGf_prime()
        self.calc_dG0_path()
        
        
    @property
    def p_h2(self):
        """Get the partial pressure of hydrogen."""
        return self._pH2

    def set_p_h2(self, value):
        """Set the partial pressure of hydrogen."""
        self._pH2 = value
    
    @property
    def p_co2(self):
        """Get the partial pressure of CO2."""
        return self._pCO2

    def set_p_co2(self, value):
        """Set the partial pressure of CO2."""
        self._pCO2 = value
        
    @property
    def T(self):
        """Get the temperature."""
        return self._T

    def set_T(self, value):
        """Set the temperature."""
        self._T = value
        #if the temperature is changed: update dGf' values and dg0 of reactions
        self.get_dGf_prime()
        self.calc_dG0_path()
        
    @property
    def maxPi(self):
        """Get the Pi pool."""
        return self._maxPi

    def set_maxPi(self, value):
        """Set the Pi pool."""
        self._maxPi = value
        
    @property
    def maxCoA(self):
        """Get the Pi pool."""
        return self._maxCoA

    def set_maxCoA(self, value):
        """Set the Pi pool."""
        self._maxCoA = value
        
    @property
    def dGatp(self):
        """Get the value of dGatp."""
        return self._dGatp

    def set_dGatp(self, value):
        """Set the value of dGatp."""
        self._dGatp = value
    
    def printreactions(self):
        """ Print all pathway reactions, to see and check if the stoichiometric matrix was set up correctly. """
        equations = {}
        
        for i, enz in enumerate(self._reactions):
            S = self._stoich[:,i]
            sub = ''
            prod = ''
            for j in range(len(S)):
                 if S[j] < 0:
                     sub += f'{abs(S[j])} {self._compounds[j]} + '
                 if S[j] > 0:
                    prod += f' {abs(S[j])} {self._compounds[j]} + '
                    
            sub = sub[0:-2]   #remove extra plus at the end for both sides of the reaction
            prod = prod[0:-2]
                
            eq = sub + '<-->' + prod
            equations[enz] = eq
        
        for key in equations.keys():
            print(f'\n {key}: {equations[key]}')
        
        self._equations = equations
        
        return self._equations 
    
    def get_dGf_prime(self):
        """ Get physiological delta G of formation for all compounds.
            Code (adapted) from eQ: https://equilibrator.readthedocs.io/en/latest/equilibrator_examples.html#Using-formation-energies-to-calculate-reaction-energies
            Coupled to database, accounted for pH and temperature."""
        
        cc = ComponentContribution()

        # changing the aqueous environment parameters
        cc.p_h = Q_(self._p_h)              
        # cc.p_mg = Q_(3.0)                 # is default
        # cc.ionic_strength = Q_("0.25M")   # is default
        cc.temperature = Q_(f"{self._T}K")    
        
        # obtain a list of compound objects using `get_compound`
        compound_list = [cc.get_compound(f"{c_id}") for c_id in list(self._compound_ids)]

        # appply standard_dg_formation on each one, and pool the results in 3 lists
        standard_dgf_mu, sigmas_fin, sigmas_inf = zip(*map(cc.standard_dg_formation, compound_list))
        standard_dgf_mu = np.array(standard_dgf_mu)
        sigmas_fin = np.array(sigmas_fin)
        sigmas_inf = np.array(sigmas_inf)

        # we now apply the Legendre transform to convert from the standard ΔGf to the standard ΔG'f
        delta_dgf_list = np.array([
            cpd.transform(cc.p_h, cc.ionic_strength, cc.temperature, cc.p_mg).m_as("kJ/mol")
            for cpd in compound_list ])
        standard_dgf_prime_mu = standard_dgf_mu + delta_dgf_list
        
        #store physiological deltaG of formation of compounds
        self._dGfprime = standard_dgf_prime_mu
        
        #now account for ratio compounds: generate 'relative' dGf'
        ratios = ['rNADH', 'rNADPH', 'rFd', 'rATP']
        #ids of oxidized versions e-carries + ADP, same order as in ratios
        accompanying_c_id = ['bigg.metabolite:nad', 'bigg.metabolite:nadp', 'kegg:C00139', 'bigg.metabolite:adp']
        
        # obtain a list of complemnentary compound objects (ratio compounds) using `get_compound`
        r_compounds_list = [cc.get_compound(f"{c_id}") for c_id in accompanying_c_id]

        # appply standard_dg_formation on each one, and pool the results in 3 lists
        dgf_mu_rc, sigmas_fin_rc, sigmas_inf_rc = zip(*map(cc.standard_dg_formation, r_compounds_list))
        dgf_mu_rc = np.array(dgf_mu_rc)
        sigmas_fin_rc = np.array(sigmas_fin_rc)
        sigmas_inf_rc = np.array(sigmas_inf_rc)

        # we now apply the Legendre transform to convert from the standard ΔGf to the standard ΔG'f
        delta_dgf_r_list = np.array([
            cpd.transform(cc.p_h, cc.ionic_strength, cc.temperature, cc.p_mg).m_as("kJ/mol")
            for cpd in r_compounds_list ])
        dgf_prime_rc = dgf_mu_rc + delta_dgf_r_list
        
        #loop through compounds and check if they are a ratio of two compounds
        for i, c in enumerate(self._compounds_copy):
            if c in ratios:
                #if so: subtract the dGf' of the accompanying compound to get a relative dGf'
                self._dGfprime[i] += - dgf_prime_rc[ratios.index(c)]
                
        return self._dGfprime
    
    def calc_dG0_path(self):
        """     Function that calculates and saves the deltaG0 values of the reactions.
                This is automatically called when a new pathway is initialized. """
                
        #first calculate dg0 values based on compound formation energies provided
        self._dg0 = self._stoich_copy.T @ self._dGfprime
        
        #save or get the dg0 of hydrogenase
        if 'hyd' in self._reactions:
            i_hyd = self._reactions.index('hyd')
            
            #normalize to production of 1 hydrogen: 2 fdred- + 2 H+ --> 2fdox + H2
            i_H2 = self._compounds.index('H2')
            self._dg0_hyd = self._dg0[i_hyd]/self._stoich[i_H2, i_hyd]
            
        else:
            self._dg0_hyd = self.get_hyd_dg0()
            
        #save all dg0 values as attribute of the pathway object
        return self._dg0, self._dg0_hyd
    
    def check_element_balance(self):
        """     Function that checks the element and charge balance of the reactions.
                This is automatically called when a new pathway is initialized. 
                Unless all balances close, you cannot proceed."""
                
        #make empty matrix (amount of elements + charge by amount of reactions) to store balance
        check = np.empty([self._element_comp.shape[1], len(self._reactions)])
        for i in range(self._element_comp.shape[1]):
            balance = self._element_comp[:,i] @ self._stoich
            check[i,:] = balance
        
        #save balance as attribute of class
        self._balance_result = check
        
        if not np.all(check == 0):
            raise ValueError(
                "One or more of the element and charge balances does not close! Check your reactions.")
        
        return self._balance_result
    

    def get_hyd_dg0(self):
        #hydrogenase reaction
        # 2 fdred- + 2 H+ --> 2fdox + H2
        #components = [Fdred, Fdox, H2]
        hyd_comps = ['kegg:C00138', 'kegg:C00139', 'bigg.metabolite:h2']
        
        hyd_S = [-2, 2, 1]
        
        cc = ComponentContribution()

        # changing the aqueous environment parameters
        cc.p_h = Q_(self._p_h)              
        # cc.p_mg = Q_(3.0)                 # is default
        # cc.ionic_strength = Q_("0.25M")   # is default
        cc.temperature = Q_(f"{self._T}K")    
        
        # obtain a list of compound objects using `get_compound`
        compound_list = [cc.get_compound(f"{c_id}") for c_id in list(hyd_comps)]

        # appply standard_dg_formation on each one, and pool the results in 3 lists
        standard_dgf_mu, sigmas_fin, sigmas_inf = zip(*map(cc.standard_dg_formation, compound_list))
        standard_dgf_mu = np.array(standard_dgf_mu)
        sigmas_fin = np.array(sigmas_fin)
        sigmas_inf = np.array(sigmas_inf)

        # we now apply the Legendre transform to convert from the standard ΔGf to the standard ΔG'f
        delta_dgf_list = np.array([
            cpd.transform(cc.p_h, cc.ionic_strength, cc.temperature, cc.p_mg).m_as("kJ/mol")
            for cpd in compound_list ])
        standard_dgf_prime_mu = standard_dgf_mu + delta_dgf_list
        
        #store physiological deltaG of formation of compounds
        hyd_dGfs = standard_dgf_prime_mu
        
        #calculate dg0' of reaction
        dG_hyd0 = hyd_S @ hyd_dGfs
        
        return dG_hyd0
