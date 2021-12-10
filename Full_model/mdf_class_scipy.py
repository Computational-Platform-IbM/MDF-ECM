# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 14:23:24 2021

@author: marit
"""
#%%
import numpy as np
from scipy.optimize import minimize

from datafile import (
    def_c_max,
    def_c_min,
    F,
    R)

from pathway_class_scipy_cc import Pathway_cc
from mdf_result_class_scipy import MDF_Result
import warnings

#%%
class MDF_Analysis(Pathway_cc):
    
    def get_constraints(self, user_defined_rNADH = False):
        
        def con_coApool(ln_conc):
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
           
            #constraint: the sum of all compounds carrying CoA should be equal to the total pool
            return totalCoA - self._maxCoA

        def con_Pipool(ln_conc):
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
            
            #constraint: the sum of all compounds carrying phosphate should be equal to the total pool (includes free Pi)
            return totalPi - self._maxPi

        def con_H2O(ln_conc):
            #get index of compound H2O
            i_H2O = self._compounds.index('H2O')
            
            #constraint: [H2O] = 1 M
            return np.exp(ln_conc[i_H2O]) - 1

        # def con_H2(ln_conc):
        #     #get index of compound H2
        #     i_H2 = self._compounds.index('H2')
            
        #     #Henry's law: p_i = H_i * c_i
        #     #with units of H in [l*atm/mol]
        #     #partial pressure proportional to mol fraction in liquid
        #     #assumption: ideal mixture, low values of x_i
        #     H_H2 = 1228         #l*atm/mol
        #     cH2 = self._pH2/H_H2
        #     #cH2 = 1e-3
            
        #     #constraint: [H2] is determined by Henry's law and hydrogen partial pressure
        #     return np.exp(ln_conc[i_H2]) - cH2
        
        # def con_CO2(ln_conc):
        #     #get index of compound CO2
        #     i_CO2 = self._compounds.index('CO2')
            
        #     #Henry's law: p_i = H_i * c_i
        #     #with units of H in [l*atm/mol]
        #     #partial pressure proportional to mol fraction in liquid
        #     #assumption: ideal mixture, low values of x_i
        #     H_CO2 = 1/0.037          #l*atm/mol
        #     cCO2 = self._pCO2/H_CO2
            
        #     ##TODO: smoothness of constraints? if/abs statements can create problems
        #     if cCO2 > def_c_max:
        #         cCO2 = def_c_max
            
        #     #constraint: [H2] is determined by Henry's law and hydrogen partial pressure
        #     return np.exp(ln_conc[i_CO2]) - cCO2

        def con_NADH_EP(ln_conc):
            i_rNADH = self._compounds.index('rNADH')
            
            #NAD+ + 2e- + H+ --> NADH
            
            n       = 2
            E0      = -320e-3
            dG0_NADH  = (-n*F*E0 )/1000  
            
            #values from Buckel & Thauer 2013
            Eprime  = -280e-3                       #V (J/C)
            
            dG_NADHprime = -n*F*Eprime/1000           #kJ
            rNADH_val = np.exp((dG_NADHprime - dG0_NADH)/(R*self._T))
            
            rNADH = ln_conc[i_rNADH]
            
            #store rNADH value from electron potential
            self._rNADH = rNADH_val
            
            return rNADH - np.log(rNADH_val)
        
        def con_NADH_set(ln_conc):
            i_rNADH = self._compounds.index('rNADH')
            rNADH = ln_conc[i_rNADH]
            
            #Get the value that was set as rNADH from the attribute rNADH of the object
            rNADH_val = self._rNADH
            
            return rNADH - np.log(rNADH_val)

        def con_NADPH(ln_conc):
            i_rNADPH = self._compounds.index('NADPH')    
            rNADPH = ln_conc[i_rNADPH]
           
            rNADPH_val = self._rNADPH
            return rNADPH - np.log(rNADPH_val)

        def con_Fd(ln_conc):
            # Relate ratio of ferredoxin directly to hydrogen production
            # Purpose of ferredoxin in cell: hydrogen production for electron sink
            # 2 Fdred- + 2H+ --> 2 Fdox + H2
            # Not too much energy available: state dG = 0
            
            #Get index of ferredoxin and value during optimization from ln_conc array
            i_rFd = self._compounds.index('rFd')
            rFd = ln_conc[i_rFd]
            
            ##TODO: show equation where this (x and y etc) is derived from
            x = np.exp(- self._dg0_hyd/(R*self._T)) * ( 1 / np.exp( ln_conc[self._compounds.index('H2')] )  )
            y = x**(1/2)
            rFd_val = 1/y
            
            #save value as attribute
            self._rFd = rFd_val
        
            return rFd - np.log(rFd_val)

        #create list of constraints
        cons = [{'type': 'eq', 'fun': con_coApool},
                {'type': 'eq', 'fun': con_Pipool},
                {'type': 'eq', 'fun': con_H2O}]
        
        #conditional additions to constraints
        #if rNADH is NOT defined by user, but based on the electron potential
        if user_defined_rNADH == False:
            cons += [{'type': 'eq', 'fun': con_NADH_EP}]
        else:
            cons += [{'type': 'eq', 'fun': con_NADH_set}]
        #if the following compounds are actually in the metabolic network:
        if 'rNADPH' in self._compounds:
            cons += [{'type': 'eq', 'fun': con_NADPH}]
        
        #Lines below commented for now: not using partial pressure! It's about the liquid phase
        # if 'H2' in self._compounds:
        #     cons += [{'type': 'eq', 'fun': con_H2}]
        # if 'CO2' in self._compounds:
        #     cons += [{'type': 'eq', 'fun': con_CO2}]
        
        if 'rFd' in self._compounds:
           cons += [{'type': 'eq', 'fun': con_Fd}]
            
        return cons
    
    def get_bounds(self, phys_bounds = False):
        """ Function to the component bounds for the MDF optimization problem. """
        
        #create bounds list
        bnds = []

        #components that are either fixed, or ratios, that don't need to abide to the default bounds of phys. concentrations
        comp_exceptions = ['H+', 'H2O', 'H2', 'rFd', 'rNADH', 'rNADPH', 'rATP']

        #loop through compounds
        for i in range(self._Nc):
            #check for exceptions
            if self._compounds[i] in comp_exceptions:
                bnds += [(None, None)]
            else:
                #take default min and max concentrations; physiological boundaries
                if phys_bounds == True:
                    bnds += [(np.log(def_c_min), np.log(def_c_max))]
                if phys_bounds == False:
                    bnds += [(None, None)]
        
        return bnds
    
    
    def execute_mdf_basis(self, set_fixed_c=False, user_defined_rNADH = False, phys_bounds = False):
        """     Function to optimize for the MDF of the pathway.
                Can be done with or without fixed concentrations for specific compounds. 
                Default is without fixed concentrations or pathway energy.   """
        
        def max_mdf(ln_conc):  
            i_Pi = self._compounds.index('Pi')
            cH = 10**-self._p_h
            
            rATP = np.exp( (self._dGatp - self._dGatp0) / (R*self._T)) * np.exp(ln_conc[i_Pi]) * (cH)
            
            #set dg_prime as a function of dg0_prime and variables ln_conc
            dg_prime = self._dg0 + ( R*self._T * self._stoich.T @ ln_conc ) + ( R * self._T * self._rATP_in_reaction * np.log(rATP) )
                
            #the optimization target is the minimum driving force (mdf) of the pathway
            #mdf is the reaction with the least amount of driving force (so the highest value)
            mdf = max(dg_prime)
            return mdf
        
        #get bounds for scipy minimize
        bnds = self.get_bounds(phys_bounds)
        
        #if you want to fix product/substrate concentrations, update bounds for that compound
        if set_fixed_c == True:
            #loop through array of fixed_c to check for specified concentrations
            for i in range(len(self._fixed_c)):
                if not np.isnan(self._fixed_c[i]):
                    #set bound as very small window of [specified value , specified value + 1e-9]
                    bnds[i] = (np.log(self._fixed_c[i]), np.log(self._fixed_c[i]+1e-9))
        
        #convert bounds list to tuple; required for scipy.optimize.minimize
        bnds = tuple(bnds)
        
        #get constraints for scipy minimize
        cons = self.get_constraints(user_defined_rNADH)
        
        #initial values
        conc0 = [np.log(1e-4)] * self._Nc
        #in E.coli, free Pi = 10 mM
        i_Pi = self._compounds.index('Pi')
        conc0[i_Pi] = np.log(10e-3)
        
        #call minimizer
        res = minimize(max_mdf, conc0, #method='SLSQP', 
                       tol=self._tol_conc, bounds = bnds, constraints = cons, options = {'maxiter': 2000})

        #check if optimizer succeeded
        if res.success == False:
            #try again with less tolerance (for example 1e-9 --> 1e-8)
            self.set_solver_tol(self._tol_conc*10)   
           
            res = minimize(max_mdf, conc0, #method='SLSQP', 
                           tol=self._tol_conc, bounds = bnds, constraints = cons, options = {'maxiter': 2000})
            
            #if optimizer did not succeed again: raise error with cause of optimization failure
            if res.success == False:
                raise ValueError(
                    f"Optimization failed: {res.message}")
        
        #get results
        opt_conc = np.exp(res.x)
        #dg_prime_opt = np.zeros(len(self._reactions))

        rATP = np.exp( (self._dGatp - self._dGatp0) / (R*self._T)) * opt_conc[i_Pi] * (10**-self._p_h)
        dg_prime_opt = self._dg0 + ( R*self._T * self._stoich.T @ res.x ) + ( R * self._T * self._rATP_in_reaction * np.log(rATP) )
        
        i_rNADH = self._compounds.index('rNADH')
        self._rNADH = opt_conc[i_rNADH]
        
        if 'rNADPH' in self._compounds:
            self._rNADPH = opt_conc[self._compounds.index('rNADPH')]
        
        #Check results: sum dgs equal to overall dg?
        sum_dg = sum(dg_prime_opt)
        
        overall_dg0 = self._S_netR_copy.T @ self._dGfprime
        overall_dg_prime = overall_dg0 + ( R*self._T * self._S_netR.T @ res.x ) + ( R * self._T * self._netATP * np.log(rATP) )
        
        #Floor values so that float precision does not matter as much
        check = np.floor(sum_dg) - np.floor(overall_dg_prime)
        #If the difference between the two floored values is not zero, something is wrong: raise error
        if check != 0:
            raise ValueError(
                "The sum of reaction energies is not equal to the overall reaction energy!")
       

        #create instance of MDF result class
        return MDF_Result(opt_conc, 
                          dg_prime_opt, 
                          overall_dg_prime,
                          self._dg0, 
                          self._reactions, 
                          self._compounds, 
                          self._S_netR, 
                          self._rATP_in_reaction,
                          self._T, 
                          self._p_h, 
                          self._pH2, 
                          self._pCO2,
                          self._maxCoA, 
                          self._maxPi, 
                          self._rNADH, 
                          self._rNADPH,
                          self._rFd,
                          self._dGatp,
                          self._dGatp0)
    
    @property
    def solver_tol(self):
        """Get the solver tolerance."""
        return self._tol_conc

    def set_solver_tol(self, value):
        """Set the solver tolerance."""
        self._tol_conc = value
    
    @property
    def rNADH(self):
        """Get the NADH/NAD+ ratio."""
        return self._rNADH

    def set_rNADH(self, value):
        """Set the NADH/NAD+ ratio."""
        self._rNADH = value
        
    @property
    def rNADPH(self):
        """Get the NADPH/NADP+ ratio."""
        return self._rNADPH

    def set_rNADPH(self, value):
        """Set the NADPH/NADP+ ratio."""
        self._rNADPH = value
    
    def to_default(self):
        """ Set all values to default values again """
        
    def fix_compound_value(self, comp: str, conc: float):
        if conc > 10e-2 or conc < 1e-6:
            warnings.warn('Warning: a fixed concentration is outside of the physiological boundaries!')
            
        i_comp = self._compounds.index(comp)
        self._fixed_c[i_comp] = conc
        
        return self._fixed_c

