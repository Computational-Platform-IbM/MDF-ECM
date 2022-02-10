# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 10:10:52 2021

@author: marit
"""

#%%
import numpy as np
import copy

from MEPfunctions import importpath 

from datafile import (
    F,
    default_T,
    default_pH,
    default_pH2,
    default_pCO2,
    dGatp0,
    default_dGatp,
    default_Pipool,
    default_CoApool,
    default_EP_NADH,
    default_EP_NADPH,
    default_EP_Fd)

import scipy.stats
from equilibrator_api import ComponentContribution, Q_
from typing import Dict, List

import warnings
#disable runtimewarnings for clarity in output
warnings.filterwarnings("ignore", category=RuntimeWarning) 

#%%

class Pathway_cc(object):
    
    def __init__(self, filename, T = default_T, pH = default_pH):
        """Create a Pathway object with default settings, importing the Pathway from an Excel sheet. """
        
        #initialize pathway with default values
        self._pH        = pH  
        self._T         = T     #K
            
        #set default values for CoA and Pi pool
        self._maxPi     = default_Pipool
        self._maxCoA    = default_CoApool
        
        #set default values for NADH/NAD+, NADPH/NADP+ and Fd_red/Fd_ox
        self._rNADH     = None
        self._rNADPH    = None
        self._rFd       = None
        
        self._EP_NADH   = default_EP_NADH
        self._EP_NADPH  = default_EP_NADPH
        self._EP_Fd     = default_EP_Fd
        
        #set default values for ATP production
        self._dGatp0    = dGatp0
        self._dGatp     = default_dGatp
        
        
        #dGprime of reactions used to determine rFd
        #hydrogenase, HydABC complex, Rnf, Nfn
        #Rnf needs to have dG that allows for the production of 1/4 ATP per proton/electron
        #Rnf is normalized to 1 NADH: 2 electrons
        self._dGprime_hyd       = 0                 #kJ/mol
        self._dGprime_hydABC    = 0                 #kJ/mol
        self._dGprime_Rnf       = -self._dGatp/2    #kJ/mol
        self._dGprime_Nfn       = 0                 #kJ/mol
        
        
        #set default tolerance for solver
        self._tol_conc    = 1e-9
        
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
         self._rel_flux)        = importpath('\\' + filename)

        
        #check for reactions that have all coefficients zero and remove them
        self.check_empty_reactions()
        
        #create copy of list with all compounds and stoich matrix
        self._compounds_copy    = self._compounds.copy()
        self._stoich_copy       = self._stoich.copy()
        
        #save compound names that have a fixed concentration in the excel sheet
        self._fixed_c_names     = [self._compounds[i] for i in range(0, len(self._compounds)) if not np.isnan(self._fixed_c[i])]
        
        #default values
        self._pH2               = default_pH2       # atm
        self._pCO2              = default_pCO2      # atm
        
        #checks if the stoichiometry closes and automatically calculates dg0 values of all reactions
        self.check_element_balance()
        self.get_dGf_prime()
        self.calc_dG0_path()
        self.get_dG_rFd()
        
        #check which reactions should be included from the optimization, as they have only fixed components
        self.exclude_from_opt()
        
        #after calculating dg0 of the pathway reactions, remove rATP from the stoichiometric matrix etc
        #save involvement of ATP/ADP in separate array
        i_rATP                  = self._compounds.index('rATP')
        self._rATP_in_reaction  = self._stoich[i_rATP,:]
        self._netATP            = self._S_netR[i_rATP]
        self._netATP_copy       = self._netATP.copy()
        self._S_netR_copy       = self._S_netR.copy()

        #remove rATP from matrix and arrays
        self._compounds.pop(i_rATP)
        self._stoich            = np.delete(self._stoich, i_rATP, axis=0)
        self._fixed_c           = np.delete(self._fixed_c, i_rATP)
        self._element_comp      = np.delete(self._element_comp, i_rATP, axis=0)
        self._S_netR            = np.delete(self._S_netR, i_rATP)
        
        self._fixed_c_copy      = list(self._fixed_c)
        self._init_values       = {'T': T, 'pH': pH, 'original fixed conc': self._fixed_c_copy}
        
        #get number of compounds and reactions in pathway
        self._Nc, self._Nr  = self._stoich.shape
        
    @property
    def p_h(self):
        """Get the pH."""
        return self._pH

    def set_ph(self, value):
        """Set the pH."""
        self._pH = value
        
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
        """Get the CoA pool."""
        return self._maxCoA

    def set_maxCoA(self, value):
        """Set the CoA pool."""
        self._maxCoA = value
        
    @property
    def dGprime_hyd(self):
        """Get the dG of hydrogenase."""
        return self._dGprime_hyd

    def set_dGprime_hyd(self, value):
        """Set the dG of hydrogenase."""
        self._dGprime_hyd = value
        
    @property
    def dGprime_hydABC(self):
        """Get the dG of hydrogenase."""
        return self._dGprime_hydABC

    def set_dGprime_hydABC(self, value):
        """Set the dG of hydrogenase."""
        self._dGprime_hydABC = value
        
    @property
    def dGatp(self):
        """Get the value of dGatp."""
        return self._dGatp

    def set_dGatp(self, value):
        """Set the value of dGatp."""
        self._dGatp = value
        
    @property
    def rNADH(self):
        """Get the NADH/NAD+ ratio."""
        return self._rNADH

    def set_rNADH(self, value):
        """Set the NADH/NAD+ ratio."""
        self._rNADH = value
    
    @property
    def EP_NADH(self):
        """Get the E' of NADH. """
        return self._EP_NADH
    
    def change_EP_NADH(self, value):
        """Change the E' of NADH. """
        self._EP_NADH = value
        
    @property
    def rNADPH(self):
        """Get the NADPH/NADP+ ratio."""
        return self._rNADPH

    def set_rNADPH(self, value):
        """Set the NADPH/NADP+ ratio."""
        self._rNADPH = value
    
    @property
    def EP_NADPH(self):
        """Get the E' of NADPH. """
        return self._EP_NADPH
    
    def change_EP_NADPH(self, value):
        """Change the E' of NADPH. """
        self._EP_NADPH = value
        
    @property
    def rFd(self):
        """Get the Fd_red-/Fd_ox ratio."""
        return self._rFd

    def set_rFd(self, value):
        """Set the Fd_red-/Fd_ox ratio."""
        self._rFd = value
        
    @property
    def EP_Fd(self):
        """Get the E' of Ferredoxin. """
        return self._EP_Fd
    
    def change_EP_Fd(self, value):
        """Change the E' of Ferredoxin. """
        self._EP_Fd = value
    
    def get_dGf_prime(self):
        """ Get physiological delta G of formation for all compounds.
            Code (adapted) from eQ: https://equilibrator.readthedocs.io/en/latest/equilibrator_examples.html#Using-formation-energies-to-calculate-reaction-energies
            Coupled to database, accounted for pH and temperature."""
        
        cc = ComponentContribution()

        # changing the aqueous environment parameters
        cc.p_h              = Q_(self._pH)              
        # cc.p_mg           = Q_(3.0)                   # is default
        # cc.ionic_strength = Q_("0.25M")               # is default
        cc.temperature      = Q_(f"{self._T}K")    
        
        # obtain a list of compound objects using `get_compound`
        compound_list = [cc.get_compound(f"{c_id}") for c_id in list(self._compound_ids)]
        
        # appply standard_dg_formation on each one, and pool the results in 3 lists
        standard_dgf_mu, sigmas_fin, sigmas_inf = zip(*map(cc.standard_dg_formation, compound_list))
        standard_dgf_mu                         = np.array(standard_dgf_mu)
        sigmas_fin                              = np.array(sigmas_fin)
        sigmas_inf                              = np.array(sigmas_inf)

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
        dgf_mu_rc                               = np.array(dgf_mu_rc)
        sigmas_fin_rc                           = np.array(sigmas_fin_rc)
        sigmas_inf_rc                           = np.array(sigmas_inf_rc)

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
        
        #save the relative formation energies
        if 'rNADH' in self._compounds:
            self._dGf_rNADH     = self._dGfprime[self._compounds_copy.index('rNADH')]
        if 'rNADPH' in self._compounds:
            self._dGf_rNADPH    = self._dGfprime[self._compounds_copy.index('rNADPH')]
        if 'rFd' in self._compounds:
            self._dGf_rFd       = self._dGfprime[self._compounds_copy.index('rFd')]
                
        return self._dGfprime
    
    def calc_dG0_path(self):
        """     Function that calculates and saves the deltaG0 values of the reactions.
                This is automatically called when a new pathway is initialized. """
                
        #first calculate dg0 values based on compound formation energies provided
        self._dg0 = self._stoich_copy.T @ self._dGfprime
        
        #save all dg0 values as attribute of the pathway object
        return self._dg0
    
    def get_dG_rFd(self):
        """     Function that determines how rFd is calculated.
                Depending on the reactions present in the pathway, rFd is coupled to one of: 
                    Rnf-complex, Nfn-complex, hydrogenase or HydABC.
                
                To this end, the dG0 values of these reactions are stored, so that they can later be used to calculate rFd."""
        
        #get index of rNADH in compounds as some reactions are normalized to rNADH
        i_rNADH = self._compounds_copy.index('rNADH')
        
        #if Rnf complex in system: use this to determine rFd
        #save the dg0 of Rnf complex if this reaction is in the matrix
        if 'Rnf' in self._reactions:
            i_Rnf = self._reactions.index('Rnf')
            
            #normalize to production of 1 NADH: 2 Fdred- + NAD+ + H+ --> 2 Fdox + NADH + pmf
            self._dg0_Rnf = self._dg0[i_Rnf]/abs(self._stoich_copy[i_rNADH, i_Rnf])
        
        #save the dg0 of hydABC complex if this reaction is in the matrix
        if 'hydABC' in self._reactions or 'HydABC' in self._reactions:
            if 'hydABC' in self._reactions:
                i_change = self._reactions.index('hydABC')
                self._reactions[i_change] = 'HydABC'
            
            i_hydABC = self._reactions.index('HydABC')
            
            #normalize to production of 1 NAD: 2 Fdred- + NADH + 3H+ --> 2 Fdox + NAD+ + 2 H2
            self._dg0_hydABC = self._dg0[i_hydABC]/abs(self._stoich_copy[i_rNADH, i_hydABC])    
            
        #save the dg0 of Nfn complex if this reaction is in the matrix
        if 'Nfn' in self._reactions:
            i_Nfn = self._reactions.index('Nfn')
            
            #normalize to production of 1 NAD: 2 Fdred- + NADH + 2 NADP+ --> 2 Fdox + NAD+ + 2 NADPH
            self._dg0_Nfn = self._dg0[i_Nfn]/abs(self._stoich_copy[i_rNADH, i_Nfn])
            
        #save or get the dg0 of hydrogenase
        if 'hyd' in self._reactions:
            i_hyd = self._reactions.index('hyd')
            
            #normalize to production of 1 hydrogen: 2 fdred- + 2 H+ --> 2fdox + H2
            i_H2 = self._compounds.index('H2')
            self._dg0_hyd = self._dg0[i_hyd]/abs(self._stoich_copy[i_H2, i_hyd])
        
        rFd_reactions = ['HydABC', 'Rnf', 'Nfn', 'hyd']
        
        #save how rFd will be determined - this attribute is used in the ferredoxin constraint from mdf_class
        if any(reaction in self._reactions for reaction in rFd_reactions):
            #rFd will be determined based on one of the reactions stated above, but will be determined (tbd) in mdf_class
            self._reaction_for_rFd = 'tbd in constraint'
        else:
            #set attribute to 'None' if these reactions are not in the pathway - rFd will then be based on E'
            self._reaction_for_rFd = None
            
        return 
    
    def check_empty_reactions(self):
        #create empty list to store indices of empty columns
        empty = []
        
        #loop through all columns (reactions) to check for reactions that have all coefficients zero
        for i in range(0, self._stoich.shape[1]):
            if np.all(self._stoich[:,i] == 0):
                empty += [i]
                
        #delete empty columns from stoich matrix, from reactions list and rel_flux array
        self._stoich        = np.delete(self._stoich, empty, axis=1)
        self._reactions     = [r for r in self._reactions if self._reactions.index(r) not in empty]
        self._rel_flux      = np.delete(self._rel_flux, empty)
        
        #now check for compounds that are not used
        #update empty list to store indices of empty rows
        empty = []
        
        #loop through all rows (compounds) to check for compounds that are not used/produced in any reaction
        for j in range(0, self._stoich.shape[0]):
            if self._compounds[j] != 'rATP' and self._compounds[j] != 'Pi':
                if np.all(self._stoich[j,:] == 0):
                    empty += [j]
        
        self._stoich        = np.delete(self._stoich, empty, axis=0)
        self._compounds     = [c for c in self._compounds if self._compounds.index(c) not in empty]
        self._element_comp  = np.delete(self._element_comp, empty, axis=0)
        self._fixed_c       = np.delete(self._fixed_c, empty) 
        self._compound_ids  = np.delete(self._compound_ids, empty) 
        self._S_netR        = np.delete(self._S_netR, empty)
                
        return 
    
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
        
        #matrix should be all (almost) equal to zero
        if not np.all(check <= 10e-12):
            raise ValueError(
                "One or more of the element and charge balances does not close! Check your reactions.")
        
        return self._balance_result
    
    def exclude_from_opt(self):
        #if a reaction purely uses conserved moieties/components that have fixed concentrations, it should be removed from the optimization
        #otherwise the optimization can't push past certain reations that are fixed: a different concentration profile could be obtained
        
        #get indices of fixed values (e-carriers and components with fixed concentrations)
        fixed = self._fixed_c_names + ['rFd', 'rNADH', 'rNADPH']
        i_fixed = [i for i, c in enumerate(self._compounds) if c in fixed]
        
        i_exclude_reactions = []
        
        #for each reaction, check if the only components participating are fixed
        #get indices of participating compounds: 
        for i in range(0, self._stoich.shape[1]):
            i_participating_compounds = [j for j, c in enumerate(self._compounds) if self._stoich[j,i] != 0]
            
            #check if all of them in i_fixed: if so, exclude from optimization
            if set(i_participating_compounds).issubset(set(i_fixed)):
                i_exclude_reactions += [i]
        
        self._excl_reactions_opt = i_exclude_reactions
        
        return

    def get_hyd_dg0(self):
        #hydrogenase reaction
        # 2 fdred- + 2 H+ --> 2fdox + H2
        #components = [Fdred, Fdox, H2]
        hyd_comps = ['kegg:C00138', 'kegg:C00139', 'bigg.metabolite:h2']
        # 2 fdred- + 2 H+ --> 2fdox + H2
        hyd_S = [-2, 2, 1]
        
        cc = ComponentContribution()

        # changing the aqueous environment parameters
        cc.p_h          = Q_(self._pH)              
        cc.temperature  = Q_(f"{self._T}K")    
        
        # obtain a list of compound objects using `get_compound`
        compound_list   = [cc.get_compound(f"{c_id}") for c_id in list(hyd_comps)]

        # appply standard_dg_formation on each one, and pool the results in 3 lists
        standard_dgf_mu, sigmas_fin, sigmas_inf = zip(*map(cc.standard_dg_formation, compound_list))
        standard_dgf_mu                         = np.array(standard_dgf_mu)
        sigmas_fin                              = np.array(sigmas_fin)
        sigmas_inf                              = np.array(sigmas_inf)

        # we now apply the Legendre transform to convert from the standard ΔGf to the standard ΔG'f
        delta_dgf_list = np.array([
            cpd.transform(cc.p_h, cc.ionic_strength, cc.temperature, cc.p_mg).m_as("kJ/mol")
            for cpd in compound_list ])
        standard_dgf_prime_mu = standard_dgf_mu + delta_dgf_list
        
        #store physiological deltaG of formation of compounds
        hyd_dGfs    = standard_dgf_prime_mu
        
        #calculate dg0' of reaction
        dG_hyd0     = hyd_S @ hyd_dGfs
        
        return dG_hyd0
    
    def to_default(self):
        """ Set all values to default/initial values again """
        #Check if temperature and pH have been changed
        #if so: adjust them and update dGf' values
        #updating dGf' is conditional because this takes some time, undesirable to do this everytime!
        if self._T != self._init_values['T'] or self._pH != self._init_values['pH']:
            self._T = self._init_values['T']
            self._pH = self._init_values['pH']
            
            self.get_dGf_prime()
            self.calc_dG0_path()
        
        #ensure fixed concentrations are back to the way they are in the excel sheet
        self._fixed_c = np.array(self._init_values['original fixed conc'])
        
        #set these parameters back to default values
        self._dGatp     = default_dGatp
        self._maxPi     = default_Pipool
        self._maxCoA    = default_CoApool
        
        self._netATP = self._netATP_copy
        
        self._dGprime_hyd       = 0                 #kJ/mol
        self._dGprime_hydABC    = 0                 #kJ/mol
        self._dGprime_Rnf       = -self._dGatp/2    #kJ/mol
        self._dGprime_Nfn       = 0                 #kJ/mol
        
        self._tol_conc          = 1e-9
        ##TODO: electron carriers?
        
        return
        
    def fix_compound_value(self, comp: str, conc: float):
        if conc > 10e-2 or conc < 1e-6:
            warnings.warn('Warning: a fixed concentration is outside of the physiological boundaries!')
            
        i_comp = self._compounds.index(comp)
        self._fixed_c[i_comp] = conc
        
        return self._fixed_c
