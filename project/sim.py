from constants import *
from xdot import xdot

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt



# maximum step size for ode solver
tstep = 5*Tg/1e3

# time intervals
tspan0 = [0, 5*Tg]
tspan1 = [5*Tg, 10*Tg]
tspan2 = [10*Tg, 15*Tg]
tspan3 = [15*Tg, 20*Tg]

# initial values
# thetag, id, iq, zd, zq, vcd, vcq, igd, igq
x0 = 0, 1.41, 0, 1.25, 0, 1.96, -0.15, 1.4, -.05

# simulate!
results0 = solve_ivp(lambda t, x: xdot(t, x, V, kp, ki, R, B0, R0, wg, V_base, XL_base, w_base, 1e6/S_base, 0/S_base), tspan0, x0, max_step=tstep)
y0 = results0.y

x1 = y0[:,y0.shape[1]-1]
results1 = solve_ivp(lambda t, x: xdot(t, x, V, kp, ki, R, B0, R0, wg, V_base, XL_base, w_base, 0.5e6/S_base, 0/S_base), tspan1, x1, max_step=tstep)
y1 = results1.y

x2 = y1[:,y1.shape[1]-1]
results2 = solve_ivp(lambda t, x: xdot(t, x, V, kp, ki, R, B0, R0, wg, V_base, XL_base, w_base, 0.5e6/S_base, 200e3/S_base), tspan2, x2, max_step=tstep)
y2 = results2.y

x3 = y2[:,y2.shape[1]-1]
results3 = solve_ivp(lambda t, x: xdot(t, x, V, kp, ki, R, B0, R0, wg, V_base, XL_base, w_base, 0.25e6/S_base, 200e3/S_base), tspan3, x3, max_step=tstep)
y3 = results3.y

# concatenate results
t = np.concatenate((results0.t, results1.t, results2.t, results3.t))
results = list()
for i in range(0,len(y0)):
    data = np.concatenate((y0[i], y1[i], y2[i], y3[i]))
    results.append(data)

# plot results
for i, angle in enumerate(results[0]):
    results[0][i] = np.mod(angle, 2*pi)

P = list()
Q = list()

for i, cur in enumerate(results[0]):
    P.append((V / V_base) * results[7][i])
    Q.append((V / V_base) * results[8][i])

plt.figure(1)
plt.plot(t, results[0], color='b')
plt.xlabel('Time (sec)')
plt.ylabel('Angle mod 2*pi (rad)')
plt.title('Grid Angle')

plt.figure(2)
plt.plot(t, results[1], label='id',color='b')
plt.plot(t, results[2], label='iq',color='r')
plt.xlabel('Time (sec)')
plt.ylabel('dq-Current (pu)')
plt.legend(('id','iq'))
plt.title("Per Unit Currents")

plt.figure(3)
plt.plot(t, P, label='P', color='b')
plt.plot(t, Q, label='Q', color='r')
plt.xlabel('Time (sec)')
plt.ylabel('Power (pu)')
plt.legend(('P','Q'))
plt.title("Real and Reactive Power")

plt.figure(4)
plt.plot(t, results[5], label='v_c_d',color='b')
plt.plot(t, results[6], label='v_c_q',color='r')
plt.xlabel('Time (sec)')
plt.ylabel('Voltage (pu)')
plt.legend(('v_c_d','v_c_q'))
plt.title("Per Unit LCL Capacitor Voltage")

plt.figure(5)
plt.plot(t, results[7], label='i_g_d',color='b')
plt.plot(t, results[8], label='i_g_q',color='r')
plt.xlabel('Time (sec)')
plt.ylabel('dq-Current (pu)')
plt.legend(('i_g_d','i_g_q'))
plt.title("Per Unit Grid Currents")

# for i in range(0, 9):
#     plt.figure(i)
#     plt.plot(t,results[i])
plt.show()