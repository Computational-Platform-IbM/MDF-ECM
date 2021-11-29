# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 10:10:52 2021

@author: marit
"""

#%%
import numpy as np
import sys

#if MEPfunctions.py not in same folder, add the required folder to path
sys.path.append('C:\\Users\marit\Documents\LST\MSc\MEP\Scipy MDF\MDF-ECM')

from MEPfunctions import importpath 

from datafile import (
    F,
    default_T,
    default_pH,
    default_pH2)


#%%

class Pathway(object):
    
    def __init__(self, filename):
        """Create a Pathway object with default settings, importing the Pathway from an Excel sheet. """
        
        #initialize pathway with default values
        self._p_h   = default_pH  
        self._T     = default_T     #K
        self._pH2   = default_pH2   #atm
        
        #set default values for CoA and Pi pool
        self._maxPi     = 20e-3     #M
        self._maxCoA    = 10e-3     #M
        
        #set default values for NADH/NAD+ and NADPH/NADP+
        self._rNADH     = 0.05
        self._rNADPH    = 100
        
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
         self._deltaGf0, 
         self._S_netR, 
         self._stoich, 
         self._rel_flux)     = importpath('\\' + filename)
        
        #get number of compounds and reactions in pathway
        self._Nc, self._Nr  = self._stoich.shape
        
        #checks if the stoichiometry closes and automatically calculates dg0 values of all reactions
        self.check_element_balance()
        self.calc_dG0_path()
        
        
    @property
    def p_h(self):
        """Get the pH."""
        return self._p_h

    def set_p_h(self, value):
        """Set the pH."""
        self._p_h = value
        
    @property
    def p_h2(self):
        """Get the partial pressure of hydrogen/concentration."""
        return self._pH2

    def set_p_h2(self, value):
        """Set the partial pressure of hydrogen/concentration."""
        self._pH2 = value
        
    @property
    def T(self):
        """Get the temperature."""
        return self._T

    def set_T(self, value):
        """Set the temperature."""
        self._T = value
        
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
    
    def calc_dG0_path(self):
        """     Function that calculates and saves the deltaG0 values of the reactions.
                This is automatically called when a new pathway is initialized. """
                
        #first calculate dg0 values based on compound formation energies provided
        self._dg0 = self._stoich.T @ self._deltaGf0

        #then correct for electron carriers involved + H2 production
        #deltag0 determined based on electron potentials

        #loop through reactions
        for i in range(self._stoich.shape[1]):
            #check if NADH/NAD+ participates in matrix
            if 'rNADH' in self._compounds:
                #values from Buckel & Thauer 2013
                # NAD+ + 2e + H+ --> NADH
                i_rNADH = self._compounds.index('rNADH')
                n       = 2
                E0      = -320e-3
                dG0_NADH  = (-n*F*E0 )/1000                 #kJ

                if self._stoich[i_rNADH, i] != 0:
                    S = self._stoich[i_rNADH, i]
                    self._dg0[i] += S*dG0_NADH
            
            #check if ferredoxin in matrix
            if 'rFd' in self._compounds:
                #values from Buckel & Thauer 2013
                # Fd_ox + 2e- --> Fd_red-2
                i_rFd = self._compounds.index('rFd')
                n       = 2
                E0      = -400e-3
                dG0_Fd  = (-n*F*E0 )/1000                   #kJ
                
                #check if Fd_red/Fd_ox participates in reaction
                if self._stoich[i_rFd, i] != 0:
                    S = self._stoich[i_rFd, i]
                    self._dg0[i] += S*dG0_Fd
            
            #check if H2 participates in matrix
            if 'H2' in self._compounds:
                #values from Buckel & Thauer 2013
                # 2H+ + 2e- --> H2
                i_H2 = self._compounds.index('H2')
                n = 2
                E0 = -414e-3
                dG0_H2 = (-n*F*E0 )/1000
                
                #check if H2 participates in reaction
                if self._stoich[i_H2, i] != 0:
                    S = self._stoich[i_H2, i]
                    self._dg0[i] += S*dG0_H2
            
        #save all dg0 values as attribute of the pathway object
        return self._dg0
    
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
    


