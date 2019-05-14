'''
Created on May 14, 2019

@author: baumj

EE 534 - HW3 - Problem 2
Write averaged converter model 
'''

import numpy as np
from scipy.integrate import solve_ivp #ODE45 
import matplotlib.pyplot as plt


#Scientific Constants
Kb = 1.38064852e-23 #Boltzman's constant J K^-1
qe = 1.6021766208e-19 #electron charge

#Constants given in problem
ig_stc = 4.6147
Rs = 0.5371
eta = 0.87223
i0 = 7.1712e-13
Rsh = 419.78
Nc = 128
Tk = 308 #temp in Kelvin (from BJ paper)

#Ns = ? #number of series modules
#Np = ? #number of parallel modules 
err = 1e-8

vt = Nc*Kb/qe*Tk*eta


def PVModel(t,x):  #x is an array of state variables [id, iq, Zd, Zq, Vdc^2, Zdc, theta^g, Zpll]
    
    
    dxdt = np.array("derivs of x") #the derivatives of the state variables
    return dxdt
