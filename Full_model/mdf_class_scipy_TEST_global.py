# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 14:23:24 2021

@author: marit
"""
#%%
import numpy as np
from scipy.optimize import differential_evolution, NonlinearConstraint

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

        def all_cons(ln_conc):
            #create empty return array for all applicable constraints, lb&ub
            cons    = []
  
            #get indices of compounds carrying phosphate
            i_p = np.where(self._element_comp[:,-2] != 0)[0]
            #if there are no compounds carrying phosphate in the matrix (left), then length will be either 0 or 1 (Pi will be there)
            #only account for phosphate pool constraint if there are compounds carrying phosphate
            if len(i_p) > 1:
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
                con_Pipool = totalPi - self._maxPi
                cons += [con_Pipool]
                
            #similar for Pi-pool: constraint for CoA-pool to be included or not?
            #get indices of compounds carrying CoA
            i_coA = []
            for i, string in enumerate(self._compounds): 
                if 'CoA' in string: 
                    i_coA += [i]
            
            if len(i_coA) > 0:
                totalCoA = 0
                #add all concentrations together
                for i in i_coA:
                    compconc = np.exp(ln_conc[i])
                    totalCoA += compconc
               
                #constraint: the sum of all compounds carrying CoA should be equal to the total pool
                con_CoApool = totalCoA - self._maxCoA
                cons += [con_CoApool]
            
            #if NADH is in the metabolic network:
            if 'rNADH' in self._compounds:
                i_rNADH = self._compounds.index('rNADH')
                rNADH = ln_conc[i_rNADH]
                
                #if rNADH is NOT defined by user, but based on the electron potential
                if user_defined_rNADH == False:
                    #NAD+ + 2e- + H+ --> NADH           
                    n       = 2

                    dG_NADHprime = -n*F*self._EP_NADH/1000           #kJ
                    rNADH_val = np.exp((dG_NADHprime - self._dGf_rNADH)/(R*self._T))
                    
                    rNADH = ln_conc[i_rNADH]
                    
                    #store rNADH value from electron potential
                    self._rNADH = rNADH_val
                    
                    con_NADH_EP = rNADH - np.log(rNADH_val)
                    #cons        += [con_NADH_EP]
                else:
                    #value used in optimization should be the same as the value of the attribute that is set by the user
                    if self._rNADH == None:
                        raise ValueError(
                            "No set value for rNADH was found.")
                    con_NADH_set = rNADH - np.log(self._rNADH)
                    #cons        += [con_NADH_set]
                    
            #if NADPH is in the metabolic network:
            if 'rNADPH' in self._compounds:
                i_rNADPH = self._compounds.index('rNADPH')    
                rNADPH = ln_conc[i_rNADPH]
                
                #when rNADPH is NOT defined by user, it's based on the electron potential ('EP')
                if user_defined_rNADPH == False:
                    #NADP+ + 2e- + H+ --> NADPH
                    n       = 2

                    dG_NADPHprime   = -n*F*self._EP_NADPH/1000           #kJ
                    rNADPH_val      = np.exp((dG_NADPHprime - self._dGf_rNADPH)/(R*self._T))
                    
                    rNADPH = ln_conc[i_rNADPH]
                    
                    #store rNADH value from electron potential
                    self._rNADPH = rNADPH_val
                    
                    con_NADPH_EP    =  rNADPH - np.log(rNADPH_val)
                    #cons            += [con_NADPH_EP]
                else:
                    #value used in optimization should be the same as the value of the attribute that is set by the user
                    if self._rNADPH == None:
                        raise ValueError(
                            "No set value for rNADPH was found.")
                    con_NADPH_set   =  rNADPH - np.log(self._rNADPH)
                    #cons            += [con_NADPH_set]
                    
            #if ferredoxin is in the metabolic network:  
            if 'rFd' in self._compounds:
                i_rFd = self._compounds.index('rFd')
                rFd = ln_conc[i_rFd]
                
                if user_defined_rFd == False:
                    if self._reaction_for_rFd:
                        rFd_dict = {}
                        ##TODO: show equation where this (x and y etc) is derived from    
                        if 'Rnf' in self._reactions:
                            #2 Fdred- + NAD+ + H+ --> 2 Fdox + NADH + pmf
                            # dG' = dG0 + R T ln(...)
                            x = np.exp( (self._dGprime_Rnf - self._dg0_Rnf) /(R*self._T)) / ( self._rNADH)  
                            y = x**(1/2)
                            rFd_dict['Rnf'] = 1/y
                            
                        if 'HydABC' in self._reactions:
                            # Relate ratio of ferredoxin directly to hydrogen production through ebif reaction (HydABC complex)
                            # Purpose of ferredoxin in cell: hydrogen production for electron sink
                            # 2 Fdred- + NADH + 3H+ --> 2 Fdox + NAD+ + 2 H2
                            # Energy available: related to self._dGprime_hydABC
                            x = np.exp( (self._dGprime_hydABC - self._dg0_hydABC) /(R*self._T)) * ( self._rNADH  / (np.exp( ln_conc[self._compounds.index('H2')] )**2)  )
                            y = x**(1/2)
                            rFd_dict['HydABC'] = 1/y
                        
                        if 'Nfn' in self._reactions:
                            x = np.exp( (self._dGprime_Nfn - self._dg0_Nfn) /(R*self._T)) * ( self._rNADH  / (self._rNADPH**2)  )
                            y = x**(1/2)
                            rFd_dict['Nfn'] = 1/y
                        
                        if 'hyd'in self._reactions:
                            # Relate ratio of ferredoxin directly to hydrogen production through hydrogenase reaction
                            # Purpose of ferredoxin in cell: hydrogen production for electron sink
                            # 2 Fdred- + 2H+ --> 2 Fdox + H2
                            # Energy available: related to self._dGprime_hyd
                            x = np.exp( (self._dGprime_hyd - self._dg0_hyd) /(R*self._T)) * ( 1 / np.exp( ln_conc[self._compounds.index('H2')] )  )
                            y = x**(1/2)
                            rFd_dict['hyd'] = 1/y
                        
                        #determine rFd:
                        #all reactions described above require reduced ferredoxin, so..
                        #the reaction that requires the most reduction power to make it feasible, determines rFd
                        #so the max value in the dictionary
                        rFd_val = max(rFd_dict.values())
                        #save the name of the reaction determining rFd
                        self._reaction_for_rFd = list(rFd_dict.keys())[list(rFd_dict.values()).index(rFd_val)]

                    #if there is no reaction to determine rFd (see pathway.get_dG_rFd())
                    #so if self._reaction_for_rFd = None
                    else:
                        # Ratio based on electron potential
                        # Fd_ox + + e- --> Fd_red-
                        n       = 1
                        
                        dG_Fdprime = -n*F*self._EP_Fd/1000           #kJ
                        rFd_val = np.exp((dG_Fdprime - self._dGf_rFd)/(R*self._T))
                    
                    #save value as attribute
                    self._rFd = rFd_val
                    
                else:
                    #Get the value that was set as rFd from the attribute rFd of the object
                    rFd_val = self._rFd
                    
                con_rFd  =   rFd - np.log(rFd_val)
                #cons        += [con_rFd]
                
            return cons 
        
        nlc     = NonlinearConstraint(all_cons, 0, 0)
        return nlc
    
    
    def get_bounds(self, phys_bounds = False):
        """ Function to the component bounds for the MDF optimization problem. """
        
        #create bounds list
        bnds = []

        #components that are either fixed, or ratios, that don't need to abide to the default bounds of phys. concentrations
        comp_exceptions = ['H2O', 'rFd', 'rNADH', 'rNADPH']

        #loop through compounds
        for i in range(self._Nc):
            #check for exceptions
            if self._compounds[i] in comp_exceptions:
                if self._compounds[i] == 'H2O':
                    bnds.append((np.log(1), np.log(1)))
                else:
                    bnds.append((np.log(1e-4), np.log(1e-5)))
            else:
                #take default min and max concentrations; physiological boundaries
                if phys_bounds == True:
                    bnds.append((np.log(def_c_min), np.log(def_c_max)))
                if phys_bounds == False:
                    bnds.append((np.log(1e-12), 10000000))
            
        return bnds
    
    
    def execute_mdf_basis(self, set_fixed_c=False, user_defined_rNADH = False, user_defined_rNADPH = False, user_defined_rFd = False, phys_bounds = False):
        """     Function to optimize for the MDF of the pathway.
                Can be done with or without fixed concentrations for specific compounds. 
                Default is without fixed concentrations or pathway energy.   """
        
        
        #function that is minimization target
        def max_mdf(ln_conc):  
            i_Pi = self._compounds.index('Pi')
            cH = 10**-self._pH
            
            rATP = np.exp( (self._dGatp - self._dGatp0) / (R*self._T)) * np.exp(ln_conc[i_Pi]) * (cH)
            
            #set dg_prime as a function of dg0_prime and variables ln_conc
            dg_prime = self._dg0 + ( R*self._T * self._stoich.T @ ln_conc ) + ( R * self._T * self._rATP_in_reaction * np.log(rATP) )
            
            print('conc')
            print(np.exp(ln_conc))
            print('dg')
            print(dg_prime)
            #the optimization target is the minimum driving force (mdf) of the pathway
            #mdf is the reaction with the least amount of driving force (so the highest value)
            #don't look at the reactions that are fixed, so ignore indices that are in self._excl_reactions_opt
            mdf_check = [val for i, val in enumerate(dg_prime) if i not in self._excl_reactions_opt]
            #check maximum value: that is the minimization target
            mdf = max(mdf_check)
            
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
        
        #get constraints for scipy minimize
        cons = self.get_constraints(user_defined_rNADH, user_defined_rNADPH, user_defined_rFd)
        
        #c0 = np.array([np.log(1e-4)] * self._Nc)
        #call minimizer
        res = differential_evolution(max_mdf, bounds = bnds, constraints = cons, disp=True, maxiter=10)#, tol = 0.0001, x0=c0)
        
        #check if optimizer succeeded
        if res.success == False:
            #try again with less tolerance (for example 1e-9 --> 1e-8)
            self.set_solver_tol(self._tol_conc*10)   
           
            res = differential_evolution(max_mdf, bounds = bnds, constraints = cons, disp=True, maxiter=10)#, tol = 0.0001, x0=c0)
            print(res)
            #if optimizer did not succeed again: raise error with cause of optimization failure
            if res.success == False:
                raise ValueError(
                    f"Optimization failed: {res.message}")
        
        #get results
        opt_conc = np.exp(res.x)
        
        i_Pi = self._compounds.index('Pi')
        
        #calculate the rATP ratio from the dGprime and dG0 values and [Pi] - the latter which comes from the optimization
        rATP = np.exp( (self._dGatp - self._dGatp0) / (R*self._T)) * opt_conc[i_Pi] * (10**-self._pH)
        dg_prime_opt = self._dg0 + ( R*self._T * self._stoich.T @ res.x ) + ( R * self._T * self._rATP_in_reaction * np.log(rATP) )
        
        #get e-carrier ratios from optimization and save
        if 'rNADH' in self._compounds:
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
       
        plot_dg0 = self._dg0
        
        #account for extra ATP production in system through pmf
        if 'Rnf' in self._reactions:
            i_Rnf = self._reactions.index('Rnf')
            #ATP produced from pmf = energy available in Rnf divided by the energy that it costs to generate ATP
            #add this to the value that was already there from SLP reactions
            self._ATP_pmf = dg_prime_opt[i_Rnf]/-self._dGatp
        else:
            self._ATP_pmf = None

        #create instance of MDF result class
        return MDF_Result(opt_conc, 
                          dg_prime_opt, 
                          overall_dg_prime,
                          plot_dg0,  
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
                          self._ATP_pmf,
                          self._dGprime_hyd,
                          self._dGprime_hydABC,
                          self._excl_reactions_opt)
    
    @property
    def solver_tol(self):
        """Get the solver tolerance."""
        return self._tol_conc

    def set_solver_tol(self, value):
        """Set the solver tolerance."""
        self._tol_conc = value
    
    @property
    def init_conc(self):
        """Get the initial concentration of compounds for the optimizer."""
        return self._c0
    
    def set_init_conc(self, value):
        """Get the initial concentration of compounds for the optimizer."""
        self._c0 = value
    
    def sensitivity_analysis(self, var: str, values: List[float], vary_compound_conc = False, 
                             set_fixed_c = True, user_defined_rNADH = False, user_defined_rNADPH = False, 
                             user_defined_rFd = False, phys_bounds = True):
        
        """Perform a sensitivity analysis of the pathway.
        
        'Var' has to be one of: [T, pH, Pi, CoA, dGprime_hyd, dGprime_hydABC, dGatp, rNADH, rNADPH, rFd]. Unless 'vary_compound_conc' = True, then 'var' can be any compound participating in the reaction.
        
        'Values' has to be a list of floats that are the different values of the parameter for which the pathway will be analyzed."""
        
        #list of options of variables that can be adjusted
        options = ['T', 'pH',  'Pi-pool', 'CoA-pool', 'dGprime_hyd', 'dGprime_hydABC', 'dGatp', 'rNADH', 'rNADPH', 'rFd']
        
        result_objects = []
        
        if vary_compound_conc == False:
            if var not in options:
                raise KeyError('The variable called cannot be varied through this function.')
        
            if var == 'T':
                for i, val in enumerate(values):
                    self.set_T(val)
                    result = self.execute_mdf_basis(set_fixed_c, user_defined_rNADH, user_defined_rNADPH, user_defined_rFd, phys_bounds)
                    result_objects += [result]
                    
            if var == 'pH':
                for i, val in enumerate(values):
                    self.set_ph(val)
                    result = self.execute_mdf_basis(set_fixed_c, user_defined_rNADH, user_defined_rNADPH, user_defined_rFd, phys_bounds)
                    result_objects += [result]
                
            if var == 'Pi-pool':
                for i, val in enumerate(values):
                    self.set_maxPi(val)
                    result = self.execute_mdf_basis(set_fixed_c, user_defined_rNADH, user_defined_rNADPH, user_defined_rFd, phys_bounds)
                    result_objects += [result]
                
            if var == 'CoA-pool':
                for i, val in enumerate(values):
                    self.set_maxCoA(val)
                    result = self.execute_mdf_basis(set_fixed_c, user_defined_rNADH, user_defined_rNADPH, user_defined_rFd, phys_bounds)
                    result_objects += [result]
                
            if var == 'dGprime_hyd':
                for i, val in enumerate(values):
                    self.set_dGprime_hyd(val)
                    result = self.execute_mdf_basis(set_fixed_c, user_defined_rNADH, user_defined_rNADPH, user_defined_rFd, phys_bounds)
                    result_objects += [result]
                
            if var == 'dGprime_hydABC':
                for i, val in enumerate(values):
                    self.set_dGprime_hydABC(val)
                    result = self.execute_mdf_basis(set_fixed_c, user_defined_rNADH, user_defined_rNADPH, user_defined_rFd, phys_bounds)
                    result_objects += [result]
            
            if var == 'dGatp':
                for i, val in enumerate(values):
                    self.set_dGatp(val)
                    result = self.execute_mdf_basis(set_fixed_c, user_defined_rNADH, user_defined_rNADPH, user_defined_rFd, phys_bounds)
                    result_objects += [result]
                
            if var == 'rNADH':
                for i, val in enumerate(values):
                    self.set_rNADH(val)
                    result = self.execute_mdf_basis(set_fixed_c, user_defined_rNADH, user_defined_rNADPH, user_defined_rFd, phys_bounds)
                    result_objects += [result]
            
            if var == 'rNADPH':
                for i, val in enumerate(values):
                    self.set_rNADPH(val)
                    result = self.execute_mdf_basis(set_fixed_c, user_defined_rNADH, user_defined_rNADPH, user_defined_rFd, phys_bounds)
                    result_objects += [result]
            
            if var == 'rFd':
                for i, val in enumerate(values):
                    self.set_rFd(val)
                    result = self.execute_mdf_basis(set_fixed_c, user_defined_rNADH, user_defined_rNADPH, user_defined_rFd, phys_bounds)
                    result_objects += [result]
                
        else:
            comp = var
            
            for i, val in enumerate(values):
                self.fix_compound_value(comp, val)
                result = self.execute_mdf_basis(set_fixed_c, user_defined_rNADH, user_defined_rNADPH, user_defined_rFd, phys_bounds)
                result_objects += [result]
        
        if vary_compound_conc == False:
            return MDF_Sens_Analysis_Result(result_objects)
        else:
            return MDF_Sens_Analysis_Result(result_objects, comp, values)


