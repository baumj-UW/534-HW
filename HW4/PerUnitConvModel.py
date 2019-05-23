'''
Created on May 20, 2019

@author: baumj
EE534 - HW4
Recreate converter model in per unit
'''
from sympy.physics.quantum.sho1d import omega
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

## Config from Midterm and HW2 ##
L = 100e-6 #averaged model inductance (from MT)
Rac = 0.75e-3 #averaged model resistance (from MT)
C = 19.25e-3
Hc = C/2
Kdc_p = Hc*2*np.pi*10 
Kdc_i = Hc*(2*np.pi*10)**2
Hpll = 1/Vpk
Ki_pll = Hpll*(2*np.pi*500)**2 
Kp_pll = Hpll*2*2*np.pi*500  
Kp_cc = 0.1
Ki_cc = 0.75

err = 1e-8

VT = Nc*Kb/qe*Tk*eta
k2 = Rs*Rsh/(Rs+Rsh)  

## Base values for per unit conv
Sbase = 1e6
Vnom = np.sqrt(2/3)*480
Vrms = np.sqrt(2)*Vnom
Vbase = Vnom
Inom = Sbase/(3*Vnom)
Irms = np.sqrt(2)*Inom
Ibase = Inom
omega_base = 2*np.pi*60
XL_base = L*omega_base*Inom/Vnom
Zac_base = Vnom/Inom
Zdc_base = Zac_base ## not sure if this is valid...

## Per unit values
Rac_pu = Rac/Zac_base
Ki_cc_pu = Ki_cc*Inom/(omega_base*Vnom)
Kp_cc_pu = Kp_cc*Inom/Vnom
Kdc_i_pu = Kdc_i*(2/3)*Zdc_base
Kdc_p_pu = Kdc_p*(2/3)*Zdc_base
Ki_pll_pu = Ki_pll*Vbase/omega_base 
BC_pu = omega_base*C*Vnom/Inom
Vpk_pu = Vpk/Vbase



## Simulation times ##
T = 1/60 #period 
STEP1 = 5*T
STEP2 = 10*T
STEP3 = 15*T
STEP4 = 20*T
SUB_INT = 500 ## number of subintervals for timesteps
SIM_STEPS = [STEP1,STEP2,STEP3,STEP4]

def PVconvModel(t,x,vref_dc,Ns,Np,ig,Rsh,Rs,k2,id_0,Qref,\
                Rac_pu,omega_base,XL_base,Ki_cc_pu,BC_pu,\
                Kdc_i_pu,Ki_pll_pu,Kp_cc_pu,Vpk_pu,Kdc_p_pu):  
    #x is an array of state variables [id, iq, Zd, Zq, Vdc^2, Zdc, theta^g, Zpll]
    #all state variables and internal vars in PU
    
    #i_dq = dq current grid side of inv
    #v_dq = dq voltage on grid side of inv
    #Zdq = integrator states = "Ceff voltages"
    #Zdc = integrator state in PI
    
    ## state variables
    i_d, i_q, z_d, z_q, v2dc, Zdc, thetaHat_g, Zpll = x 
    
    ## algebraic eqns ###
    v_d = Vpk_pu ## more accurate using cos(theta-theta) 
    v_q = 0 
 
    ##PV module eqns##  --> convert back to actual for this part
    Pin = PVin(np.sqrt(v2dc)*Vbase, ig, id_0, Ns, Np, Rsh, Rs, k2, i0)/Sbase
    
    # Current Controller reference
    Pref = CalcPref(Zdc, v2dc, vref_dc, Pin, Kdc_p_pu) #this calc is in pu
    iref_d, iref_q = PQac_Idq(v_d,v_q,Pref,Qref)
    
    wg_hat = Zpll +  Ki_pll_pu*v_q 
    vt_d = z_d + (iref_d - i_d)*Kp_cc_pu - XL_base*wg_hat*i_q + v_d ## control voltage input to AC side
    vt_q = z_q + (iref_q - i_q)*Kp_cc_pu + XL_base*wg_hat*i_d + v_q
    
    ## differential equations ##
#     dId_dt = (1/L)*(vt_d - i_d*Rac - v_d + L*wg_hat*i_q) 
#     dIq_dt = (1/L)*(vt_q - i_q*Rac - v_q - L*wg_hat*i_d)
    dId_dt = (omega_base/XL_base)*(vt_d - i_d*Rac_pu - v_d + XL_base*wg_hat*i_q) 
    dIq_dt = (omega_base/XL_base)*(vt_q - i_q*Rac_pu - v_q - XL_base*wg_hat*i_d) 
    
#     dZd_dt = Ki_cc*(iref_d - i_d) 
#     dZq_dt = Ki_cc*(iref_q - i_q)
    dZd_dt = omega_base*Ki_cc_pu*(iref_d - i_d) 
    dZq_dt = omega_base*Ki_cc_pu*(iref_q - i_q) 
    
#     dV2dc_dt = (2/C)*(Pin - (3/2)*(vt_d*i_d + vt_q*i_q)) 
#     dZdc_dt = Kdc_i*(v2dc - vref_dc**2)  
    dV2dc_dt = (3*omega_base/BC_pu)*(Pin - (vt_d*i_d + vt_q*i_q))  
    dZdc_dt = Kdc_i_pu*(v2dc - vref_dc**2) 
    
    #dxdt[6] = wg_hat # = Zpll +  Kpll~ *v_q # d(thetag_hat)/dt
    #dxdt[7] = Kp_pll*v_q # d(Zpll)/dt 
    
    #dxdt = [dId_dt, dIq_dt, dZd_dt, dZq_dt, dV2dc_dt, dZdc_dt, wg_hat, Ki_pll*v_q] ##the derivatives of the state eqns
    dxdt = [dId_dt, dIq_dt, dZd_dt, dZq_dt, dV2dc_dt, dZdc_dt, omega_base*wg_hat, Ki_pll_pu*v_q] ##the derivatives of the state eqns
    return dxdt

def PVin(vdc,ig,id_0,Ns,Np,Rsh,Rs,k2,i0):
    ##PV module eqns## 
    vm = vdc/Ns 
    k1 = Rsh*(vm+Rs*ig)/(Rs+Rsh)
    ipv_d, im = NewtonsMethod(g, id_0,ig,k1,k2,i0,Rsh) #returns im from NR function; 
    return vdc*im*Np 


def CalcPref(Zdc,v2dc,vref_dc,Pin,Kdc_p=Kdc_p_pu):
    return Zdc + Kdc_p*(v2dc - vref_dc**2) + Pin
    

def g(x,k1,k2,i0,vt=VT): 
    #function for calculating PV array current Iref
    return i0*(math.exp((k1-k2*x)/vt)-1) - x  

def dgdx(x,k1,k2,i0,vt=VT): 
    #derivative function for calculating PV array current Iref
    return i0*(-k2/vt)*math.exp((k1-k2*x)/vt) - 1

def NewtonsMethod(g, x, ig, k1, k2, i0, Rsh, tol=err): 
    while True:
        x1 = x - g(x,k1,k2,i0) / dgdx(x,k1,k2,i0) 
        t = abs(x1 - x)
        if t < tol:
            break
        x = x1
    im = ig-x-(k1-k2*x)/Rsh    
    return x, im

def Idq_PQac(i_d,i_q,v_d=Vpk_pu,v_q=0): #Return AC P, Q; default v_dq constant
    Vdq = [[v_d,v_q],[v_q,-v_d]]
    return np.matmul(Vdq,[[i_d],[i_q]]) ##remove 3/2 factor to return pu (according to notes)

def PQac_Idq(v_d,v_q,Pac,Qac): #Return Id, Iq from Vdq and PQ ref
    return (1)/(v_d**2 + v_q**2)*np.matmul([[v_d,v_q],[v_q,-v_d]],[[Pac],[Qac]]) ## remove (2/3) factor to return pu



## compute equilibrium  to initialize model
#init condit
id_0 = 0  

V_m = np.linspace(0,90,100) #vector of module voltages
   
####
irrad_arr = np.array([1, 1/2, 1/2, 1/4]) * ig_stc #array of insolation current for simulation 

id_arr = np.zeros((len(irrad_arr),len(V_m)))
im_arr = np.zeros((len(irrad_arr),len(V_m)))

Pm = [None]*len(irrad_arr)
Vmpp = np.zeros(irrad_arr.shape)
Impp = np.zeros(irrad_arr.shape)
Isc = np.zeros(irrad_arr.shape)
Voc = np.zeros(irrad_arr.shape)
id_mpp = np.zeros(irrad_arr.shape)
#calc IV curve and related MPPs for each ig    
for (j,ig) in enumerate(irrad_arr):    
    for (i,vm) in enumerate(V_m):
        k1 = Rsh*(vm+Rs*ig)/(Rs+Rsh)
        id_arr[j,i], im_arr[j,i] = NewtonsMethod(g,id_0,ig,k1,k2,i0,Rsh,err)
    Pm[j] = np.multiply(V_m,im_arr[j,:])
    Vmpp[j] = V_m[np.argmax(Pm[j])]
    Impp[j] = im_arr[j,np.argmax(Pm[j])]
    id_mpp[j] = id_arr[j,np.argmax(Pm[j])] ## use for init value in NR
    Isc[j] = np.max(im_arr[j,:])
    Voc[j] = V_m[np.argmin(np.abs(im_arr[j,:]))]

Ns = np.ceil(Vmpp_array/Vmpp[0]) ## choose Ns and Np based on ig=ig_stc
Np = np.round(Pmax_array/(Ns*np.max(Pm[0])))
###



Vref_dc = Vmpp*Ns ## array of ref Vdc corresponding to all simulation igs
Qref_sim = [0, 0, 200e3, 200e3] #Q input control for simulation steps

results = [None]*4
initPV = np.zeros((4,8))  # initial conditions for 4 simulation changes [id, iq, Zd, Zq, Vdc^2, Zdc, theta^g, Zpll]
initPV[0,:] = [1.6986e3/(2*Inom), 0/Ibase, 1.2740/Vbase, 0/Vbase, (Vref_dc[0]/Vbase)**2, -3.246e3/(Sbase), 0, 1]#2*np.pi*60/omega_base]
eval_times = np.array([np.linspace(0,STEP1,SUB_INT),\
                       np.linspace(STEP1,STEP2,SUB_INT),\
                       np.linspace(STEP2,STEP3,SUB_INT),\
                       np.linspace(STEP3,STEP4,SUB_INT)]) 

#Solve ODE for each eval time  ## input calculated with NR above
results[0] = solve_ivp(lambda t, x: PVconvModel(t, x, Vref_dc[0]/Vbase, Ns, Np,\
                                                irrad_arr[0], Rsh, Rs,k2,\
                                                id_mpp[0],Qref_sim[0]/Sbase, Rac_pu,\
                                                omega_base,XL_base,Ki_cc_pu,\
                                                BC_pu,Kdc_i_pu,Ki_pll_pu,\
                                                Kp_cc_pu,Vpk_pu,Kdc_p_pu),\
                    [0,STEP1],initPV[0,:],t_eval=eval_times[0])  

for i in range(1,len(results)):
    results[i] = solve_ivp(lambda t, x: PVconvModel(t, x, Vref_dc[i]/Vbase, Ns, Np,\
                                                irrad_arr[i], Rsh, Rs,k2,\
                                                id_mpp[i],Qref_sim[i]/Sbase, Rac_pu,\
                                                omega_base,XL_base,Ki_cc_pu,\
                                                BC_pu,Kdc_i_pu,Ki_pll_pu,\
                                                Kp_cc_pu,Vpk_pu,Kdc_p_pu),\
                        SIM_STEPS[i-1:i+1],results[i-1].y[:,-1],t_eval=eval_times[i])

# 
# 
# results[2] = solve_ivp(lambda t, x: PVconvModel(t, x, Vref_dc[2], Ns, Np, \
#                                                 irrad_arr[2], Rsh, Rs, k2,\
#                                                 id_mpp[2],Qref_sim[2]),\
#                     [STEP2,STEP3],results[1].y[:,-1],t_eval=eval_times[2])
# 
# 
# results[3] = solve_ivp(lambda t, x: PVconvModel(t, x, Vref_dc[3], Ns, Np, \
#                                                 irrad_arr[3], Rsh, Rs, k2,\
#                                                 id_mpp[3],Qref_sim[3]),\
#                     [STEP3,STEP4],results[2].y[:,-1],t_eval=eval_times[3])


## Solution complete! Plot all results

# Combine solution results to plot


Pin_sim = np.zeros((len(results),SUB_INT))  ## this works because all simulation step intervals are the same
Pref_sim = np.zeros((len(results),SUB_INT)) 
Pac_sim = np.zeros((len(results),SUB_INT)) 
Qac_sim = np.zeros((len(results),SUB_INT)) 
Vdc_ref_sim = np.ones((len(results),SUB_INT))
Qacref_sim = np.ones((len(results),SUB_INT))
sim_time = np.zeros((len(results),SUB_INT)) 
for (i,res) in enumerate(results):
    for j in range(SUB_INT): ##return values in per unit (res.y vals will come in pu)
        Pin_sim[i,j] = PVin(np.sqrt(res.y[4,j])*Vbase, irrad_arr[i], id_mpp[i], Ns, Np, Rsh, Rs, k2, i0)/Sbase #res.y[4,j] = Vdc^2 at step j
        Pref_sim[i,j] = CalcPref(res.y[5,j], res.y[4,j], Vref_dc[i]/Vbase, Pin_sim[i,j], Kdc_p_pu) #res.y[5,j] = Zdc
        Pac_sim[i,j], Qac_sim[i,j] = Idq_PQac(res.y[0,j], res.y[1,j]) #res.y[0:1] = i_dq 
    Vdc_ref_sim[i,:] = Vdc_ref_sim[i,:] * Vref_dc[i]/Vbase #plot in pu
    Qacref_sim[i,:] = Qacref_sim[i,:] * Qref_sim[i]/Sbase #plot in pu
    sim_time[i] = res.t


sim_time = sim_time.flatten()

Vdc_figs = plt.figure(1)
plt.plot(sim_time,Vdc_ref_sim.flatten(),label='Vdc ref',color='r')
for res in results:
    plt.plot(res.t,np.sqrt(res.y[4]),label='Vdc',color='b')
plt.grid(True)
plt.xlabel('Time (sec)')
plt.ylabel('Voltage (pu)')
plt.legend(('Vdc Ref','Vdc Actual',))
plt.title("DC-Link Voltage (Per Unit)")


idq_figs = plt.figure(2)
for res in results:
    plt.plot(res.t,res.y[0],label='id',color='b')
    plt.plot(res.t,res.y[1],label='iq',color='r')
plt.grid(True)
plt.xlabel('Time (sec)')
plt.ylabel('dq-Current (pu)')
plt.legend(('id','iq'))
plt.title("dq-Current Outputs (Per Unit)")

p_figs = plt.figure(3)
plt.plot(sim_time,Pac_sim.flatten(),label='P Actual')
plt.plot(sim_time,Pref_sim.flatten(),label='P Ref')
plt.plot(sim_time,Pin_sim.flatten(),label='P in')
plt.grid(True)
plt.xlabel('Time (sec)')
plt.ylabel('Active Power (W)')
plt.legend()
plt.title("Active Power Output")


q_figs = plt.figure(4)
plt.plot(sim_time,Qac_sim.flatten(),label='Q Actual')
plt.plot(sim_time,Qacref_sim.flatten(),label='Q Ref')
plt.grid(True)
plt.xlabel('Time (sec)')
plt.ylabel('Reactive Power (pu)')
plt.legend()
plt.title("Reactive Power Output (Per Unit)")


theta_figs = plt.figure(5)
for res in results:
    plt.plot(res.t,np.mod(res.y[6],2*np.pi),label='theta_g',color='b')
plt.grid(True)
plt.xlabel('Time (sec)')
plt.ylabel('Angle mod 2*pi (rad)')
plt.legend(('theta_g hat',))
plt.title("PLL Control Angle")

plt.show()               
