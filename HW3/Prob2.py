'''
Created on May 14, 2019

@author: baumj

EE 534 - HW3 - Problem 2
Write averaged converter model 
'''
import math
import numpy as np
from scipy.integrate import solve_ivp #ODE45 
from scipy import misc
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
Tk = 298 #Standard test conditions 25C
Vmpp_array = 1250.0 # Minimum Vmpp of PV array
Pmax_array = 1e6 #Total array power
Vpk = np.sqrt(2/3)*480 # peak line-neutral voltage in V
L = 690e-6 #averaged model inductance (from HW2)
Rac = 6e-3 #averaged model resistance (from HW2)


err = 1e-8

vt = Nc*Kb/qe*Tk*eta


def PVModel(t,x):  #x is an array of state variables [id, iq, Zd, Zq, Vdc^2, Zdc, theta^g, Zpll]
    
    #Zdq = integrator states = "Ceff voltages"
    #Zdc = integrator state in PI
    i_d, i_q, z_d, z_q, v2dc, Zdc, thetaHat_g, Zpll = x #this should work...
#     i_q = x[1]
#     z_d = x[2]
#     z_q = x[3]
    
    
    iref_d, iref_q = ## solve for iref_dq w/ newton raphson at each time step. I think these are from prob 1. convert to dq
    vt_d = z_d + (iref_d - i_d)*Kp - L*wg_hat*i_q + v_d ## not sure about this
    #, vt_q = ## find from new i_dq
    
    dId_dt = (1/L)*(vt_d - i_d*Rac - v_d + L*wg_hat*i_q) #vt_d and vd are different...?
    dIq_dt = (1/L)*(vt_q - i_q*Rac - v_q - L*wg_hat*i_d) #vt_q and vq are different...?
    
    dZd_dt = Ki_*(iref_d - i_d) # Ki from which controller?
    dZq_dt = Ki_*(iref_q - i_q) # Ki from which controller?
    
    dV2dc_dt = (2/C)*(Pin - (2/3)*(vt_d*i_d + vt_q*i_q)) 
    dZdc_dt = Kdc_*(v2dc - v2dc_ref) #v2dc ref from pv newton ??
    
    dxdt[6] = wg_hat # = Zpll +  Kpll~ *v_q # d(thetag_hat)/dt
    dxdt[7] = Kpll*v_q # d(Zpll)/dt  not sure where vq comes from...
    
    dxdt = np.array("derivs of x") #the derivatives of the state eqns
    return dxdt


#function for calculating PV array current Iref
def g(i_d):
    return i0*(math.exp((k1-k2*i_d)/vt)-1) - i_d

#x0 = newton(g, id_0, fprime=None, args=(),tol=err, fprime2=None ) #args might be constants to feed into g func

def NewtonsMethod(g, x, tol):#, i0,k1,k2,vt):    
    while True:
        x1 = x - g(x) / misc.derivative(g, x)
        t = abs(x1 - x)
        if t < tol:
            break
        x = x1
    im = ig-x-(k1-k2*x)/Rsh    
    return x, im
 

#init condit
id_0 = 0  #FIGURE OUT VARIABLE SETUP

V_m = np.linspace(0,90,100) #vector of module voltages
   
####
irrad_arr = np.linspace(0,ig_stc,10) #array of insolation current from 0<ig<ig_stc

id_arr = np.zeros((len(irrad_arr),len(V_m)))
im_arr = np.zeros((len(irrad_arr),len(V_m)))

k2 = Rs*Rsh/(Rs+Rsh)    
    
    
for (j,ig) in enumerate(irrad_arr):    
    for (i,vm) in enumerate(V_m):
        k1 = Rsh*(vm+Rs*ig)/(Rs+Rsh)
        id_arr[j,i], im_arr[j,i] = NewtonsMethod(g,id_0,err)

# def getPVratings(ig):
#     
#     return 

## calculations based on ig=ig_stc
Pm = np.multiply(V_m,im_arr[-1,:])
Vmpp = V_m[np.argmax(Pm)]
Impp = im_arr[-1,np.argmax(Pm)]
Isc = np.max(im_arr[-1,:])
Voc = V_m[np.argmin(np.abs(im_arr[-1,:]))]

Ns = np.ceil(Vmpp_array/Vmpp)
Np = np.round(Pmax_array/(Ns*np.max(Pm)))
###


### Part A ###
## compute equilibrium values for all states between 2V<vdc<=Voc
V_dc = np.linspace(2*Vpk,Voc)
V_m = V_dc/Ns



print("working?")