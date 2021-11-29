# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 10:18:23 2021

@author: marit
"""

#%%
def_c_min   = 10**-6                    #M
def_c_max   = 10**-2                    #M

#gas constant
R       = 8.31446261815324e-3           #kJ/(K * mol)
#Faraday constant
F       = 96485.322                     #C/mol e-

default_T       = 293.15                    #K
default_pH      = 7
default_pH2     = 1                 #atm
#cH      = 10**-pH                   #M

#fixed concentrations; cell biology in numbers (c_conmoi)
#rather fix ratios, but the book results in very different ratios?
#con_moi = ['ATP', 'ADP', 'Pi', 'NAD', 'NADH', 'NADP', 'NADPH', 'CoA']
con_moi = ['ATP', 'NADH','NADPH']

c_conmoi = {'ATP':      9.6e-3,   #M
            'ADP':      0.55e-3, 
            'Pi':       20e-3,          #source: Bionumbers
            'NAD':      2.6e-3,
            'NADH':     0.083e-3,
            'NADP':     0.0021e-3,
            'NADPH':    0.12e-3,
            'CoA':      10e-3}#,           #source: Noor 2014 (PLOS)
            #'Fd_ox':    float('NaN'),
            #'Fd_red':   float('NaN')} 
            