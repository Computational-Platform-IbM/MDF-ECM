# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 12:51:00 2021

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
class Equilibrium_Analysis(Pathway_cc):
    
    def find_conc(ln_conc):
        
        return