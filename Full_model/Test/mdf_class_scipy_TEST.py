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
        
        def con_H2O(ln_conc):
            #get index of compound H2O
            i_H2O = self._compounds.index('H2O')
            
            #constraint: [H2O] = 1 M
            return np.exp(ln_conc[i_H2O]) - 1
        
        def con_NADH_EP(ln_conc):
            i_rNADH = self._compounds.index('rNADH')
            
            #NAD+ + 2e- + H+ --> NADH           
            n       = 2

            dG_NADHprime = -n*F*self._EP_NADH/1000           #kJ
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

            dG_NADPHprime   = -n*F*self._EP_NADPH/1000           #kJ
            rNADPH_val      = np.exp((dG_NADPHprime - self._dGf_rNADPH)/(R*self._T))
            
            rNADPH = ln_conc[i_rNADPH]
            
            #store rNADH value from electron potential
            self._rNADPH = rNADPH_val
            
            return rNADPH - np.log(rNADPH_val)

        def con_NADPH_set(ln_conc):
            i_rNADPH = self._compounds.index('rNADPH')    
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
        
            return rFd - np.log(rFd_val)
        
        def con_Fd_set(ln_conc):
            i_rFd = self._compounds.index('rFd')
            rFd = ln_conc[i_rFd]
            
            #Get the value that was set as rFd from the attribute rFd of the object
            rFd_val = self._rFd
            
            return rFd - np.log(rFd_val)

        #create list of constraints, based on what compounds are part of the pathway
        #start with empty constraint list
        cons = []
        
        if 'H2O' in self._compounds:
            cons += [{'type': 'eq', 'fun': con_H2O}]
        
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
          
        # cons += [{'type': 'eq', 'fun': con_Pi}, 
        #          {'type': 'eq', 'fun': con_CoA}]
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
            elif self._compounds[i] == 'Pi':
                bnds += [(np.log(self._cPi),    np.log(self._cPi+1e-12))]
            elif self._compounds[i] == 'CoA':
                bnds += [(np.log(self._cCoA),   np.log(self._cCoA+1e-12))]
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
        
        totalCoA    = 10
        totalPi     = 10
        
        self._cPi   = 1
        self._cCoA  = 1
        
        while totalCoA > self._maxCoA or totalPi > self._maxPi:
            
            #function that is minimization target
            def max_mdf(ln_conc):  
                cH = 10**-self._pH
                rATP = np.exp( (self._dGatp - self._dGatp0) / (R*self._T)) * self._cPi * (cH)
                
                #normalise reactions in stoichiometric matrix all to 1 substrate
                norm = np.amax(self._stoich, axis=0)
                S = self._stoich/norm
                rATP_in_r = self._rATP_in_reaction/norm
                
                #set dg_prime as a function of dg0_prime and variables ln_conc
                dg_prime = self._dg0 + ( R*self._T * S.T @ ln_conc ) + ( R * self._T * rATP_in_r * np.log(rATP) )
                    
                #the optimization target is the minimum driving force (mdf) of the pathway
                #mdf is the reaction with the least amount of driving force (so the highest value)
                #don't look at the reactions that are fixed, so ignore indices that are in self._excl_reactions_opt
                mdf_check = [val for i, val in enumerate(dg_prime) if i not in self._excl_reactions_opt]
                #check maximum value: that is the minimization target
                mdf = max(mdf_check)
                
                #mintarget = abs(mdf) 
                return mdf
                #return mintarget
            
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
            if type(self._c0) == float:
                conc0 = [np.log(self._c0)] * self._Nc
            elif len(self._c0) != self._Nc:
                raise ValueError('Amount of values in c0 does not match the amount of compounds in the pathway.')
            else:
                conc0 = [np.log(self._c0)]
                
            #call minimizer
            res = minimize(max_mdf, conc0, #method='SLSQP', 
                           tol=self._tol_conc, bounds = bnds, constraints = cons, options = {'disp': True, 'maxiter': 2000})
    
            #check if optimizer succeeded
            if res.success == False:
                #try again with more tolerance (for example 1e-9 --> 1e-8)
                self.set_solver_tol(self._tol_conc*10)   
               
                res = minimize(max_mdf, conc0, #method='SLSQP',
                               tol=self._tol_conc, bounds = bnds, constraints = cons, options = {'disp': True, 'maxiter': 2000})
                
                #if optimizer did not succeed again: raise error with cause of optimization failure
                if res.success == False:
                    raise ValueError(
                        f"Optimization failed: {res.message}")
            
            #get results
            opt_conc = np.exp(res.x)
            
            norm = np.amax(self._stoich, axis=0)
            S = self._stoich/norm
            
            rATP_in_r = self._rATP_in_reaction/norm
            
            #calculate the rATP ratio from the dGprime and dG0 values and [Pi] - the latter which comes from the optimization
            rATP = np.exp( (self._dGatp - self._dGatp0) / (R*self._T)) * self._cPi * (10**-self._pH)
            dg_prime_opt = self._dg0 + ( R*self._T * S.T @ res.x ) + ( R * self._T * rATP_in_r * np.log(rATP) )
            
            dg_prime_opt = dg_prime_opt * norm
            
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
            
            # #Floor values so that float precision does not matter as much
            # check = np.floor(sum_dg) - np.floor(overall_dg_prime)
            # #If the difference between the two floored values is not zero, something is wrong: raise error
            # if check != 0:
            #     raise ValueError(
            #         "The sum of reaction energies is not equal to the overall reaction energy!")
           
            plot_dg0 = self._dg0
            
            #account for extra ATP production in system through pmf
            if 'Rnf' in self._reactions:
                i_Rnf = self._reactions.index('Rnf')
                #ATP produced from pmf = energy available in Rnf divided by the energy that it costs to generate ATP
                #add this to the value that was already there from SLP reactions
                self._ATP_pmf = dg_prime_opt[i_Rnf]/-self._dGatp
            else:
                self._ATP_pmf = None
                
            #get indices of compounds carrying CoA
            i_coA = []
            for i, string in enumerate(self._compounds): 
                if 'CoA' in string: 
                    i_coA += [i]
                    
            totalCoA = 0
            #add all concentrations together
            for i in i_coA:
                compconc = np.exp(res.x[i])
                totalCoA += compconc
            
            
            #get indices of compounds carrying phosphate
            i_p = np.where(self._element_comp[:,-2] != 0)[0]
            
            totalPi = 0
            for i in i_p:
                #don't include ATP in phosphate constraint
                if not self._compounds[i] == 'rATP]':
                    #x_p is the amount of phosphate a compound carries
                    x_p = self._element_comp[i,-2]
                    #multiply concentration by amount of phosphate in molecule
                    compconc = x_p * np.exp(res.x[i])
                    #add moles of phosphate together
                    totalPi += compconc
            
            if totalCoA > self._maxCoA:
                newCoA  = self._cCoA - self._cCoA*((totalCoA - self._maxCoA)/totalCoA)
                
                if not str(newCoA)[:14] == str(self._cCoA)[:14]:
                    self._cCoA = newCoA
                else:
                    self._cCoA -= 1e-9
            if totalPi > self._maxPi:
                newPi  = self._cPi - self._cPi*((totalPi - self._maxPi)/totalPi)
                
                if not str(newPi)[:14] == str(self._cPi)[:14]:
                    self._cPi = newPi
                else:
                    self._cPi -= 1e-9
        
        print(f'{(totalPi/self._maxPi)*100:.2f} % of Pi-pool used')
        print(f'{(totalCoA/self._maxCoA)*100:.2f} % of CoA-pool used')
            
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


