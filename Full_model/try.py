# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 11:34:32 2021

@author: marit
"""
#%%
import sys
#if MEPfunctions.py and datafile.py not in same folder, add the required folder to path
sys.path.append('C:\\Users\marit\Documents\LST\MSc\MEP\Scipy MDF\MDF-ECM')

from pathway_class_scipy_cc import Pathway_cc


#%%
glu_to_ac = Pathway_cc('try')

#%%
for i, c in enumerate(glu_to_ac._compounds):
    print(f'{c}: {glu_to_ac._dGfprime[i]:.2f}')
    
print('\n\n')
for i, r in enumerate(glu_to_ac._reactions):
    print(f'{r}: {glu_to_ac._dg0[i]:.2f}')