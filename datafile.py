# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 10:18:23 2021

@author: marit
"""

#%%
def_c_min   = 1e-6                    #M
def_c_max   = 1e-2                    #M

#gas constant
R       = 8.31446261815324e-3         #kJ/(K * mol)
#Faraday constant
F       = 96485.322                   #C/mol e-

default_T       = 298.15              #K
default_pH      = 7
default_pH2     = 0.01                #atm
default_pCO2    = 0.01                #atm

dGatp0          = -8.15               #kJ/mol
default_dGatp   = 50                  #kJ/mol

default_Pipool  = 20e-3               #M
default_CoApool = 10e-3               #M

default_c0      = 1e-4                #M

#default E' values for e-carries
#values from Buckel&Thauer 2013
default_EP_NADH     = -280e-3                       #V (J/C)
default_EP_NADPH    = -380e-3                       #V (J/C)
default_EP_Fd       = -500e-3                       #V (J/C)