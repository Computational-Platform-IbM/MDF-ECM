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
from mdf_sens_analysis_result_class import MDF_Sens_Analysis_Result
import warnings
from typing import Dict, List

#%%
class MDF_Analysis(Pathway_cc):
    
    def get_constraints(self, user_defined_rNADH, user_defined_rNADPH, user_defined_rFd):
        
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
            #values from Buckel & Thauer 2013
            Eprime  = -280e-3                       #V (J/C)
            
            dG_NADHprime = -n*F*Eprime/1000           #kJ
            rNADH_val = np.exp((dG_NADHprime - self._dGf_rNADH)/(R*self._T))
            
            rNADH = ln_conc[i_rNADH]
            
            #store rNADH value from electron potential
            self._rNADH = rNADH_val
            
            return rNADH - np.log(rNADH_val)
        
        def con_NADH_set(ln_conc):
            i_rNADH = self._compounds.index('rNADH')
            rNADH = ln_conc[i_rNADH]
            
            #value used in optimization should be the same as the value of the attribute that is set by the user
            if self._rNADH == None:
                raise ValueError(
                    "No set value for rNADH was found.")
            return rNADH - np.log(self._rNADH)
        
        def con_NADPH_EP(ln_conc):
            i_rNADPH = self._compounds.index('rNADPH')
            
            #NADP+ + 2e- + H+ --> NADPH
            n       = 2
            #values from Buckel & Thauer 2013
            Eprime  = -380e-3                       #V (J/C)
            
            dG_NADPHprime   = -n*F*Eprime/1000           #kJ
            rNADPH_val      = np.exp((dG_NADPHprime - self._dGf_rNADPH)/(R*self._T))
            
            rNADPH = ln_conc[i_rNADPH]
            
            #store rNADH value from electron potential
            self._rNADPH = rNADPH_val
            
            return rNADPH - np.log(rNADPH_val)

        def con_NADPH_set(ln_conc):
            i_rNADPH = self._compounds.index('NADPH')    
            rNADPH = ln_conc[i_rNADPH]
           
            #value used in optimization should be the same as the value of the attribute that is set by the user
            if self._rNADPH == None:
                raise ValueError(
                    "No set value for rNADPH was found.")
            return rNADPH - np.log(self._rNADPH)

        def con_Fd(ln_conc):            
            #Get index of ferredoxin and value during optimization from ln_conc array
            i_rFd = self._compounds.index('rFd')
            rFd = ln_conc[i_rFd]
            
            ##TODO: show equation where this (x and y etc) is derived from           
            if 'HydABC' in self._reactions or 'hydABC' in self._reactions:
                # Relate ratio of ferredoxin directly to hydrogen production through ebif reaction (HydABC complex)
                # Purpose of ferredoxin in cell: hydrogen production for electron sink
                # 2 Fdred- + NADH + 3H+ --> 2 Fdox + NAD+ + 2 H2
                # Energy available: related to self._dGprime_hydABC
                
                #NAD+ + 2e- + H+ --> NADH           
                n       = 2
                #values from Buckel & Thauer 2013
                Eprime  = -280e-3                       #V (J/C)
                
                dG_NADHprime = -n*F*Eprime/1000           #kJ
                rNADH_val = np.exp((dG_NADHprime - self._dGf_rNADH)/(R*self._T))
                
                x = np.exp( (self._dGprime_hydABC - self._dg0_hydABC) /(R*self._T)) * ( rNADH_val  / (np.exp( ln_conc[self._compounds.index('H2')] )**2)  )
                y = x**(1/2)
                rFd_val = 1/y
                
            elif 'H2'in self._compounds:
                # Relate ratio of ferredoxin directly to hydrogen production through hydrogenase reaction
                # Purpose of ferredoxin in cell: hydrogen production for electron sink
                # 2 Fdred- + 2H+ --> 2 Fdox + H2
                # Energy available: related to self._dGprime_hyd
                x = np.exp( (self._dGprime_hyd - self._dg0_hyd) /(R*self._T)) * ( 1 / np.exp( ln_conc[self._compounds.index('H2')] )  )
                y = x**(1/2)
                rFd_val = 1/y
            
            else:
                # Ratio based on electron potential
                # Fd_ox + + e- --> Fd_red-
                n       = 1
            
                # #values from Buckel & Thauer 2013
                Eprime  = -500e-3                       #V (J/C)
                
                dG_Fdprime = -n*F*Eprime/1000           #kJ
                rFd_val = np.exp((dG_Fdprime - self._dGf_rFd)/(R*self._T))
            
            # print(rFd_val)
            #rFd_val = 1
            #save value as attribute
            self._rFd = rFd_val
        
            return rFd - np.log(rFd_val)
        
        def con_Fd_set(ln_conc):
            i_rFd = self._compounds.index('rFd')
            rFd = ln_conc[i_rFd]
            
            #Get the value that was set as rNADH from the attribute rNADH of the object
            rFd_val = self._rFd
            
            return rFd - np.log(rFd_val)

        #create list of constraints
        cons = [{'type': 'eq', 'fun': con_coApool},
                {'type': 'eq', 'fun': con_Pipool},
                {'type': 'eq', 'fun': con_H2O}]
        
        #conditional additions to constraints
        if 'rNADH' in self._compounds:
            #if rNADH is NOT defined by user, but based on the electron potential
            if user_defined_rNADH == False:
                cons += [{'type': 'eq', 'fun': con_NADH_EP}]
            else:
                cons += [{'type': 'eq', 'fun': con_NADH_set}]
        
        #if NADPH is actually in the metabolic network:
        if 'rNADPH' in self._compounds:
            #when rNADPH is NOT defined by user, it's based on the electron potential ('EP')
            if user_defined_rNADPH == False:
                cons += [{'type': 'eq', 'fun': con_NADPH_EP}]
            else:
                cons += [{'type': 'eq', 'fun': con_NADPH_set}]
        
        #constraint for rFd conditions: is Fd involved and is the rFd value manually set or no
        if 'rFd' in self._compounds:
            if user_defined_rFd == False:
                cons += [{'type': 'eq', 'fun': con_Fd}]
            else:
                cons += [{'type': 'eq', 'fun': con_Fd_set}]
                
        #TODO: Lines below commented for now: not using partial pressure! It's about the liquid phase
        # if 'H2' in self._compounds:
        #     cons += [{'type': 'eq', 'fun': con_H2}]
        # if 'CO2' in self._compounds:
        #     cons += [{'type': 'eq', 'fun': con_CO2}]
        
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
    
    
    def execute_mdf_basis(self, set_fixed_c=False, user_defined_rNADH = False, user_defined_rNADPH = False, user_defined_rFd = False, phys_bounds = False):
        """     Function to optimize for the MDF of the pathway.
                Can be done with or without fixed concentrations for specific compounds. 
                Default is without fixed concentrations or pathway energy.   """
     
        def max_mdf(ln_conc):  
            i_Pi = self._compounds.index('Pi')
            cH = 10**-self._pH
            
            rATP = np.exp( (self._dGatp - self._dGatp0) / (R*self._T)) * np.exp(ln_conc[i_Pi]) * (cH)
            
            #set dg_prime as a function of dg0_prime and variables ln_conc
            dg_prime = self._dg0 + ( R*self._T * self._stoich.T @ ln_conc ) + ( R * self._T * self._rATP_in_reaction * np.log(rATP) )
                
            #the optimization target is the minimum driving force (mdf) of the pathway
            #mdf is the reaction with the least amount of driving force (so the highest value)
            mdf = max(dg_prime)
            #print(mdf)
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
        cons = self.get_constraints(user_defined_rNADH, user_defined_rNADPH, user_defined_rFd)
        
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

        rATP = np.exp( (self._dGatp - self._dGatp0) / (R*self._T)) * opt_conc[i_Pi] * (10**-self._pH)
        dg_prime_opt = self._dg0 + ( R*self._T * self._stoich.T @ res.x ) + ( R * self._T * self._rATP_in_reaction * np.log(rATP) )
        
        if 'rNADH' in self._compounds:
            i_rNADH = self._compounds.index('rNADH')
            self._rNADH = opt_conc[i_rNADH]
        
        if 'rNADPH' in self._compounds:
            self._rNADPH = opt_conc[self._compounds.index('rNADPH')]
        
        #Check results: sum dgs equal to overall dg?
        sum_dg = sum(dg_prime_opt)
        
        overall_dg0 = self._S_netR_copy.T @ self._dGfprime
        #multiply the overall reaction energy by the moles of substrate as the net reaction is normalized to 1
        # and the sum of individual reaction energies is not!
        overall_dg_prime = (overall_dg0 + ( R*self._T * self._S_netR.T @ res.x ) + ( R * self._T * self._netATP * np.log(rATP) ) ) * abs(self._stoich[0,0])
        
        #Floor values so that float precision does not matter as much
        check = np.floor(sum_dg) - np.floor(overall_dg_prime)
        #If the difference between the two floored values is not zero, something is wrong: raise error
        if check != 0:
            raise ValueError(
                "The sum of reaction energies is not equal to the overall reaction energy!")
       
        ##TODO: disucss approach!
        #normalize all reaction energies to 1 mol of substrate - to compare scenarios
        dg_prime_opt = dg_prime_opt/abs(self._stoich[0,0])
        plot_dg0 = self._dg0/abs(self._stoich[0,0])

        #create instance of MDF result class
        return MDF_Result(opt_conc, 
                          dg_prime_opt, 
                          overall_dg_prime/abs(self._stoich[0,0]),
                          plot_dg0, #self._dg0, 
                          self._reactions, 
                          self._compounds, 
                          self._fixed_c_names,
                          self._S_netR, 
                          self._rATP_in_reaction,
                          self._T, 
                          self._pH, 
                          self._pH2, 
                          self._pCO2,
                          self._maxCoA, 
                          self._maxPi, 
                          self._rNADH, 
                          self._rNADPH,
                          self._rFd,
                          self._dGatp,
                          self._dGatp0,
                          self._netATP,
                          self._dGprime_hyd,
                          self._dGprime_hydABC)
    
    @property
    def solver_tol(self):
        """Get the solver tolerance."""
        return self._tol_conc

    def set_solver_tol(self, value):
        """Set the solver tolerance."""
        self._tol_conc = value
    
    def sensitivity_analysis(self, var: str, values: List[float], vary_compound_conc = False):
        """Perform a sensitivity analysis of the pathway.
        
        'Var' has to be one of: [T, pH, Pi, CoA, dGprime_hyd, dGprime_hydABC, dGatp, rNADH, rNADPH, rFd]. Unless 'vary_compound_conc' = True, then 'var' can be any compound participating in the reaction.
        
        'Values' has to be a list of floats that are the different values of the parameter for which the pathway will be analyzed."""
        
        #begin by setting the Pathway back to default settings every new sensitivity analysis
        self.to_default()
        
        #list of options of variables that can be adjusted
        options = ['T', 'pH',  'Pi-pool', 'CoA-pool', 'dGprime_hyd', 'dGprime_hydABC', 'dGatp', 'rNADH', 'rNADPH', 'rFd']
        
        result_objects = []
        
        if vary_compound_conc == False:
            if var not in options:
                raise KeyError('The variable called cannot be varied through this function.')
        
            if var == 'T':
                for i, val in enumerate(values):
                    self.set_T(val)
                    result = self.execute_mdf_basis(set_fixed_c=True, phys_bounds = True)
                    result_objects += [result]
                    
            if var == 'pH':
                for i, val in enumerate(values):
                    self.set_ph(val)
                    result = self.execute_mdf_basis(set_fixed_c=True, phys_bounds = True)
                    result_objects += [result]
                
            if var == 'Pi-pool':
                for i, val in enumerate(values):
                    self.set_maxPi(val)
                    result = self.execute_mdf_basis(set_fixed_c=True, phys_bounds = True)
                    result_objects += [result]
                
            if var == 'CoA-pool':
                for i, val in enumerate(values):
                    self.set_maxCoA(val)
                    result = self.execute_mdf_basis(set_fixed_c=True, phys_bounds = True)
                    result_objects += [result]
                
            if var == 'dGprime_hyd':
                for i, val in enumerate(values):
                    self.set_dGprime_hyd(val)
                    result = self.execute_mdf_basis(set_fixed_c=True, phys_bounds = True)
                    result_objects += [result]
                
            if var == 'dGprime_hydABC':
                for i, val in enumerate(values):
                    self.set_dGprime_hydABC(val)
                    result = self.execute_mdf_basis(set_fixed_c=True, phys_bounds = True)
                    result_objects += [result]
            
            if var == 'dGatp':
                for i, val in enumerate(values):
                    self.set_dGatp(val)
                    result = self.execute_mdf_basis(set_fixed_c=True, phys_bounds = True)
                    result_objects += [result]
                
            if var == 'rNADH':
                for i, val in enumerate(values):
                    self.set_rNADH(val)
                    result = self.execute_mdf_basis(set_fixed_c=True, phys_bounds = True,  user_defined_rNADH = True)
                    result_objects += [result]
            
            if var == 'rNADPH':
                for i, val in enumerate(values):
                    self.set_rNADPH(val)
                    result = self.execute_mdf_basis(set_fixed_c=True, phys_bounds = True, user_defined_rNADPH = True)
                    result_objects += [result]
            
            if var == 'rFd':
                for i, val in enumerate(values):
                    self.set_rFd(val)
                    result = self.execute_mdf_basis(set_fixed_c=True, phys_bounds = True, user_defined_rFd = True)
                    result_objects += [result]
                
        else:
            comp = var
            
            for i, val in enumerate(values):
                self.fix_compound_value(comp, val)
                result = self.execute_mdf_basis(set_fixed_c=True, phys_bounds = True)
                result_objects += [result]
        
        if vary_compound_conc == False:
            return MDF_Sens_Analysis_Result(result_objects)
        else:
            return MDF_Sens_Analysis_Result(result_objects, comp, values)


