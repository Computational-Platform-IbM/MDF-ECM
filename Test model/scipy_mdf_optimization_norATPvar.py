# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 17:17:42 2021

@author: marit
"""

#%%
from scipy.optimize import minimize
import numpy as np
#import pandas as pd
from MEPfunctions import importpath
import matplotlib.pyplot as plt

from datafile import (
    def_c_max,
    def_c_min,
    R,
    default_T,
    default_pH)

#%% Data
c_min   = def_c_min     #M
c_max   = def_c_max     #M

T       = default_T     #K
pH      = default_pH
cH      = 10**-pH       #M
p_H2    = 0.01           #atm
p_CO2   = 0.01          #atm

F       = 96485.322                     #C/mol e-

maxPi   = 20e-3
maxCoA  = 10e-3

#%% Optimization target for scipy.optimize.minimize

def max_mdf(ln_conc, dg0, rATP_in_reaction):  
    dg_prime = np.zeros(len(reactions))
    i_Pi = compounds.index('Pi')

    dGatp0  = -8.15    #kJ/mol
    dGatp   = 50
    
    rATP = np.exp( (dGatp - dGatp0) / (R*T)) * np.exp(ln_conc[i_Pi]) * (cH)
    dg_prime = dg0 + ( R*T*stoich.T @ ln_conc ) + ( R * T * rATP_in_reaction * np.log(rATP) )
    
    #the optimization target is the lowest driving force of the pathway: this should be increased (so reaction energy lowered)
    mdf = max(dg_prime)
    
    return mdf

#import pathway
(reactions, compounds, element_balance, 
 fixed_c, deltaGf0, S_netR, stoich, rel_flux) = importpath('\\glu_to_acetate_ratios_new.xlsx') #\\glu_to_ac_e-bfc_ratios_new.xlsx')  #


#%% Calculate dg0 values of reactions
#first calculate dg0 values based on compound formation energies provided
dg0 = stoich.T @ deltaGf0

#then correct for electron carriers involved + H2 production
#loop through reactions
for i in range(stoich.shape[1]):
    #check if NADH/NAD+ participates in reaction
    if 'rNADH' in compounds:
        #values from Buckel & Thauer 2013
        # NAD+ + 2e + H+ --> NADH
        i_rNADH = compounds.index('rNADH')
        n       = 2
        E0      = -320e-3
        dG0_NADH  = (-n*F*E0 )/1000                 #kJ

        if stoich[i_rNADH, i] != 0:
            S = stoich[i_rNADH, i]
            dg0[i] += S*dG0_NADH
    
    #check if ferredoxin in reaction
    if 'rFd' in compounds:
        #values from Buckel & Thauer 2013
        # Fd_ox + 2e- --> Fd_red-2
        i_rFd = compounds.index('rFd')
        n       = 2
        E0      = -400e-3
        dG0_Fd  = (-n*F*E0 )/1000                   #kJ
        
        #check if Fd_red/Fd_ox participates in reaction
        if stoich[i_rFd, i] != 0:
            S = stoich[i_rFd, i]
            dg0[i] += S*dG0_Fd
    
    #check if H2 participates in reaction
    if 'H2' in compounds:
        #values from Buckel & Thauer 2013
        # 2H+ + 2e- --> H2
        i_H2 = compounds.index('H2')
        n = 2
        E0 = -414e-3
        dG0_H2 = (-n*F*E0 )/1000
        
        #check if H2 participates in reaction
        if stoich[i_H2, i] != 0:
            S = stoich[i_H2, i]
            dg0[i] += S*dG0_H2


i_rATP = compounds.index('rATP')
#save involvement of rATP in array
rATP_in_reaction = stoich[i_rATP,:]

#remove rATP from matrix and arrays
compounds.pop(i_rATP)
stoich = np.delete(stoich, i_rATP, axis=0)
fixed_c = np.delete(fixed_c, i_rATP)
element_balance = np.delete(element_balance, i_rATP, axis=0)

#%% Bounds    
#get number of compounds and reactions in pathway
Nc, Nr = stoich.shape

#create bounds list
bnds = []

#components that are either fixed, or ratios, that don't need to abide to the default bounds of phys. concentrations
comp_exceptions = ['H+', 'H2O', 'H2', 'rFd', 'rNADH', 'rNADPH', 'rATP']

#loop through compounds
for i in range(Nc):
    #check for exceptions
    if compounds[i] in comp_exceptions:# or compounds[i] == 'Pi':
        bnds += [(None, None)]
    else:
        #take default min and max concentrations; physiological boundaries
        bnds += [(np.log(def_c_min), np.log(def_c_max))]

#%% Constraints for scipy.optimize.minimize

def con_coApool(ln_conc):
    #CoA pool constraint: totalCoA = maxCoA
    #totalCoA= total amount of compounds carrying CoA
    #maxCoA = CoA pool
    
    #maxCoA = 10e-3   #M
    i_coA = []
    for i, string in enumerate(compounds): 
        if 'CoA' in string: 
            i_coA += [i]
            
    totalCoA = 0
    for i in i_coA:
        compconc = np.exp(ln_conc[i])
        totalCoA += compconc
        
    return totalCoA - maxCoA

def con_Pipool(ln_conc):
    #phosphate pool constraint: totalPi = maxPi
    #totalPi = total amount of compounds carrying phosphate
    #maxPi = phosphate pool
    
    #maxPi = 30e-3   #M
    i_p = np.where(element_balance[:,-2] != 0)[0]
            
    totalPi = 0
    for i in i_p:
        x_p = element_balance[i,-2]
        compconc = x_p * np.exp(ln_conc[i])
        totalPi += compconc
        
    return totalPi - maxPi

def con_H(ln_conc):
    i_H = compounds.index('H+')
    return np.exp(ln_conc[i_H]) - cH

def con_H2O(ln_conc):
    i_H2O = compounds.index('H2O')
    return np.exp(ln_conc[i_H2O]) - 1

def con_H2(ln_conc):
    i_H2 = compounds.index('H2')
    
    #Henry's law: p_i = H_i * c_i
    #with units of H in [l*atm/mol]
    #partial pressure proportional to mol fraction in liquid
    #assumption: ideal mixture, low values of x_i
    H_H2 = 1228         #l*atm/mol
    cH2 = p_H2/H_H2

    return np.exp(ln_conc[i_H2]) - cH2

def con_CO2(ln_conc):
    i_CO2 = compounds.index('CO2')
    
    #Henry's law: p_i = H_i * c_i
    #with units of H in [l*atm/mol]
    #partial pressure proportional to mol fraction in liquid
    #assumption: ideal mixture, low values of x_i
    H_CO2 = 1/0.037         #l*atm/mol
    cCO2 = p_CO2/H_CO2
    
    if cCO2 > def_c_max:
        cCO2 = def_c_max

    return np.exp(ln_conc[i_CO2]) - cCO2


def con_NADH(ln_conc):
    i_rNADH = compounds.index('rNADH')
    rNADH = ln_conc[i_rNADH]

    ##TODO: option to adjust rNADH? leave free to optimize?
    rNADH_val = 0.05
    return rNADH - np.log(rNADH_val)

def con_NADPH(ln_conc):
    i_rNADPH = compounds.index('rNADPH')    
    rNADPH = ln_conc[i_rNADPH]
   
    ##TODO: option to adjust rNADPH?
    rNADPH_val = 100
    return rNADPH - np.log(rNADPH_val)

def con_Fd(ln_conc):
    i_rFd = compounds.index('rFd')
    # dG' = - n F E'
    # combine with dG' = dGo + RTln(P/S)
    # Fd_ox + 2e- --> Fd_red-2
    
    #values from Buckel & Thauer 2013
    n       = 2
    Eprime  = -500e-3                       #V (J/C)
    
    dG_Fdprime = -n*F*Eprime/1000           #kJ
    rFd_val = np.exp((dG_Fdprime - dG0_Fd)/(R*T))
    
    rFd = ln_conc[i_rFd]
    return rFd - np.log(rFd_val)

#create list of constraints
cons = [{'type': 'eq', 'fun': con_coApool},
        {'type': 'eq', 'fun': con_Pipool},
        {'type': 'eq', 'fun': con_H},
        {'type': 'eq', 'fun': con_H2O},
        {'type': 'eq', 'fun': con_NADH}]            #comment this line to optimize rNADH, instead of fixing it      

if 'rNADPH' in compounds:
    cons += [{'type': 'eq', 'fun': con_NADPH}]
if 'H2' in compounds:
    cons += [{'type': 'eq', 'fun': con_H2}]
if 'CO2' in compounds:
    cons += [{'type': 'eq', 'fun': con_CO2}]
if 'rFd' in compounds:
    cons += [{'type': 'eq', 'fun': con_Fd}]

#if you want to fix product/substrate concentrations, update bounds for that compound
for i in range(len(fixed_c)):
    if not np.isnan(fixed_c[i]):
        bnds[i] = (np.log(fixed_c[i]), np.log(fixed_c[i]+1e-9))
          
#tolerance and initial values
tol_conc = 1e-8
conc0 = [np.log(1e-4)] * Nc
conc0[-1] = 7.9
conc0[-2] = -3
#convert bounds list to tuple; required for scipy.optimize.minimize
bnds = tuple(bnds)

#call minimizer
res = minimize(max_mdf, conc0, args=(dg0, rATP_in_reaction), method='SLSQP', 
               tol=tol_conc, bounds = bnds, constraints = cons, options = {'maxiter': 2000})

#%% Get results
if res.success == False:
    raise ValueError(
        f"Optimization failed: {res.message}")
    
conc_values = np.exp(res.x)

#get optimized dg values of reactions
dGatp0  = -8.15    #kJ/mol
dGatp   = 50

i_Pi = compounds.index('Pi')
dg_prime_opt = np.zeros(len(reactions))

rATP = np.exp( (dGatp - dGatp0) / (R*T)) * conc_values[i_Pi] * (10**-pH)

for i, x in enumerate(rATP_in_reaction):
    dg_prime_opt[i] = dg0[i] + ( R*T * stoich[:,i] @ res.x) + (R * T * np.log(rATP**x))

#%% Plot results
ATP = S_netR[i_rATP]
S_netR = np.delete(S_netR, i_rATP)

sub = ''
prod = ''
for j in range(len(S_netR)):
     if S_netR[j] < 0:
         sub += f'{abs(S_netR[j]):.2f} {compounds[j]} + '
     if S_netR[j] > 0:
        prod += f' {abs(S_netR[j]):.2f} {compounds[j]} + '

if ATP > 0:
    prod += f'{ATP} ATP +'
elif ATP < 0:
    sub += f'{abs(ATP)} ATP +'
            
sub = sub[0:-2]   #remove extra plus at the end for both sides of the reaction
prod = prod[0:-2]


    
netreaction_eq = sub + '<-->' + prod

conditions = ''
for comp in comp_exceptions:
    if comp in compounds and comp != 'H2O':
        if comp == 'H+':
            conditions += f'\t pH = {pH}'.expandtabs()
        elif comp == 'H2':
            conditions += f'\t pH2 = {p_H2} atm'.expandtabs()
        else:
            i = compounds.index(comp)
            conditions += f'\t {comp} = {conc_values[i]:.3e}'.expandtabs()

conditions += f'\t pCO2 = {p_CO2} atm \t Pi-pool = {maxPi:.2e} M \t CoA-pool = {maxCoA:.2e} M'.expandtabs()



#%% check coA and Pi pool constraint
i_coA = []
for i, string in enumerate(compounds): 
    if 'CoA' in string: 
        i_coA += [i]

check_totalCoA = 0
for i in i_coA:
    compconc = conc_values[i]
    check_totalCoA += compconc
    
i_p = np.where(element_balance[:,-2] != 0)[0]
x_p = element_balance[i_p,-2]
check_totalPi = sum(x_p*conc_values[i_p])

#%% PLOT alternative
fig = plt.figure(figsize=(13,9))
plt.suptitle(netreaction_eq + '\n\n' + conditions)

#cumulative deltaG: course of pathway
sumdg0 = [0.0] + np.cumsum(dg0).tolist()
sumdg_opt = [0.0] + np.cumsum(dg_prime_opt).tolist()


#inidividual deltaG values of reactions
ax2 = fig.add_subplot(211)
ax2.title.set_text('Reaction energies')
ax2.plot(reactions, dg0, 'rs', label='dG0 [kJ/mol]', alpha=0.4)
ax2.plot(reactions, dg_prime_opt, 'bo', label = "dG' optimized [kJ/mol]")
ax2.grid()
ax2.set_ylabel('$\Delta$G [kJ/mol]')
ax2.set_xlabel('Reactions')
ax2.set_xticklabels(reactions, fontsize=8, rotation=90)
ax2.legend()

#concentrations
remove = []
for i, comp in enumerate(compounds):
    if comp in comp_exceptions:
        remove += [i]

plot_compounds = np.delete(np.array(compounds), remove)
plot_conc = np.delete(np.array(conc_values), remove)
        
ax3 = fig.add_subplot(212)
ax3.title.set_text('Compound concentrations')
ax3.plot(plot_compounds, plot_conc, 'bo')
ax3.plot(plot_compounds, plot_conc, 'b--', alpha =0.4)
ax3.grid()
ax3.set_ylabel('Concentration [M]')
ax3.set_xlabel('Compounds')
ax3.set_xticks(plot_compounds)
ax3.set_xticklabels(plot_compounds, fontsize=8, rotation=90)
ax3.set_ylim(10**-7, 10**-1)
ax3.plot([0, len(plot_conc)], [10**-6, 10**-6], 'k-')
ax3.plot([0, len(plot_conc)], [10**-2, 10**-2], 'k-')
ax3.set_yscale('log')

plt.tight_layout()

#%% PLOT
