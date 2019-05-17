'''
Created on May 12, 2019

@author: baumj

HW3 - Probem 1
Implement Newton Raphson iteration for single diode PV model
'''

import math     # Math functions 
import cmath    # Complex math function conj, rect
import numpy as np  # Methods for linear algebra
from numpy.linalg import inv
from scipy.optimize import newton
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

err = 1e-8

vt = Nc*Kb/qe*Tk*eta

#init condit
id_0 = 0  #FIGURE OUT VARIABLE SETUP

vm_arr = np.linspace(0,100,500) #vector of module voltages

def g(i_d):
    return i0*(math.exp((k1-k2*i_d)/vt)-1) - i_d


def NewtonsMethod(g, x, tol):    
    while True:
        x1 = x - g(x) / misc.derivative(g, x)
        t = abs(x1 - x)
        if t < tol:
            break
        x = x1
    im = ig-x-(k1-k2*x)/Rsh    
    return x, im

id_sol = np.zeros((vm_arr.shape))
im_sol = np.zeros((vm_arr.shape))

ig = ig_stc ## for this prob
k2 = Rs*Rsh/(Rs+Rsh)    
    
for (i,vm) in enumerate(vm_arr):
    k1 = Rsh*(vm+Rs*ig)/(Rs+Rsh)
    id_sol[i], im_sol[i] = NewtonsMethod(g,id_0,err)
    
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
fig1 = plt.figure(1) 
plt.plot(vm_arr,im_sol)    
plt.title("PV Module I-V Curve")
plt.xlabel("Voltage (V)")
plt.ylabel("Current (A)")
plt.ylim(bottom=0)
plt.xlim(left=0)
plt.grid(b=None, which='major', axis='both')
plt.annotate("Isc=%.2f" % Isc,(0,Isc),xytext=(10,4),arrowprops={'arrowstyle':'->'})
plt.annotate("Voc=%.1f" % Voc,(Voc,0),xytext=(70,0.5),arrowprops={'arrowstyle':'->'})
plt.annotate("MPP=%d W  Vmpp=%.1f V  Impp=%.2f A" % (np.max(Pm),Vmpp,Impp),(Vmpp,Impp),xytext=(40,5),arrowprops={'arrowstyle':'->'})



## Part D ## plot for whole array
irrad_arr = np.array([1, 3/4, 1/2, 0.1]) * ig_stc

id_arr = np.zeros((len(irrad_arr),len(vm_arr)))
im_arr = np.zeros((len(irrad_arr),len(vm_arr)))

k2 = Rs*Rsh/(Rs+Rsh)    
    
fig2 = plt.figure(2)
   
for (j,ig) in enumerate(irrad_arr):    
    for (i,vm) in enumerate(vm_arr):
        k1 = Rsh*(vm+Rs*ig)/(Rs+Rsh)
        id_arr[j,i], im_arr[j,i] = NewtonsMethod(g,id_0,err)

    text = "ig: " + str(np.round(ig,3))
    plt.plot(Ns*vm_arr,Np*im_arr[j,:],label=text)

   
plt.title("PV Array I-V Curves")
plt.xlabel("Voltage (V)")
plt.ylabel("Current (A)")
plt.ylim(bottom=0)
plt.xlim(left=0)
plt.grid(b=None, which='major', axis='both')
plt.legend()

plt.show()
        
