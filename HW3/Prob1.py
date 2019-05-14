'''
Created on May 12, 2019

@author: baumj

Implement Newton Raphson iteration for single diode PV model
'''

import math     # Math functions 
import cmath    # Complex math function conj, rect
import numpy as np  # Methods for linear algebra
from numpy.linalg import inv
from scipy.optimize import newton
from scipy import misc
#from scipy.integrate import odeint  #refs odeint directly instead of long pointer
#from scipy.integrate import solve_ivp #ODE45 equiv (use this instead of odeint)
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

    #Kb = 1.38064852e-23 #Boltzman's constant J K^-1
    #qe = 1.6021766208e-19 #electron charge
    
    #ig_stc = 4.6147
   # Rs = 0.5371
   # eta = 0.87223
   # i0 = 7.1712e-13
   # Rsh = 419.78
   # Nc = 128
    
    #vt = Nc*Kb/qe*Tk*eta
    #vm = Vdc/Ns  #check which version is correct





#init condit
id_0 = 0  #FIGURE OUT VARIABLE SETUP

vm_arr = np.linspace(0,88,500) #vector of module voltages

def g(i_d):#,i0,k1,k2,vt):
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

# def im(ig,id,k1,k2,Rsh):
#     return ig-id-id*(k1-k2)/Rsh

id_sol = np.zeros((vm_arr.shape))
im_sol = np.zeros((vm_arr.shape))

ig = ig_stc ## for this prob
k2 = Rs*Rsh/(Rs+Rsh)    
    
    

for (i,vm) in enumerate(vm_arr):
    k1 = Rsh*(vm+Rs*ig)/(Rs+Rsh)
    id_sol[i], im_sol[i] = NewtonsMethod(g,id_0,err)
    #x0 = newton(g, id_0, fprime=None, args=(),tol=err, fprime2=None ) #args might be constants to feed into g func
    #print(vm,sol)
    
    #print("numpy",x0)
Pm = np.multiply(vm_arr,im_sol)
Vmpp = vm_arr[np.argmax(Pm)]
Impp = im_sol[np.argmax(Pm)]
Isc = np.max(im_sol)
Voc = vm_arr[np.argmin(np.abs(im_sol))]

#label these points on the plot
print("Isc:",Isc,"A")
print("Voc:",Voc,"V")
print("MPP:",np.max(Pm),"W")
print("Vmpp:",Vmpp)
print("Impp:",Impp)

#Part C
Ns = np.ceil(1250.0/Vmpp)
Np = np.round(1e6/(Ns*np.max(Pm)))
print("Number of series modules:",Ns)
print("Number of parallel modules:",Np)
print("Total modules in array:",Np*Ns)

#plot for Part B    
plt.plot(vm_arr,im_sol)    
plt.title("PV Module I-V Curve")
plt.xlabel("Voltage (V)")
plt.ylabel("Current (A)")
plt.ylim(bottom=0)
plt.xlim(left=0)
plt.grid(b=None, which='major', axis='both')
plt.annotate("Isc=%d" % Isc,(0,Isc),xytext=(10,4),arrowprops={'arrowstyle':'->'})
plt.annotate("Voc=%d" % Voc,(Voc,0),xytext=(70,0.5),arrowprops={'arrowstyle':'->'})
plt.annotate("MPP (Vmpp,Impp)",(Vmpp,Impp),xytext=(50,4),arrowprops={'arrowstyle':'->'})
plt.show()


## Part D
irrad_arr = np.array([1, 3/4, 1/2, 0.1]) * ig_stc

id_arr = np.zeros((len(irrad_arr),len(vm_arr)))
im_arr = np.zeros((len(irrad_arr),len(vm_arr)))

k2 = Rs*Rsh/(Rs+Rsh)    
    
    
for (j,ig) in enumerate(irrad_arr):    
    for (i,vm) in enumerate(vm_arr):
        k1 = Rsh*(vm+Rs*ig)/(Rs+Rsh)
        id_arr[j,i], im_arr[j,i] = NewtonsMethod(g,id_0,err)
        #x0 = newton(g, id_0, fprime=None, args=(),tol=err, fprime2=None ) #args might be constants to feed into g func
        #print(vm,sol)
    
    text = "ig: " + str(np.round(ig,3))
    plt.plot(vm_arr,im_arr[j,:],label=text)

   
plt.title("PV Module I-V Curves")
plt.xlabel("Voltage (V)")
plt.ylabel("Current (A)")
plt.ylim(bottom=0)
plt.xlim(left=0)
plt.grid(b=None, which='major', axis='both')
plt.legend()
plt.show()
        
print("working?")
