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
# theta_g, i_d, i_q, z_d, z_q, v_c_d, v_c_q
x0 = 0, 1.6986e3/I_base, 0, 1.2740/V_base, 0/V_base, 0, 0

# simulate!
results0 = solve_ivp(lambda t, x: xdot(t, x, V, kp, ki, R, B0, wg, V_base, XL_base, w_base, 1e6/S_base, 0), tspan0, x0, max_step=tstep)
y0 = results0.y

x1 = y0[:,len(y0)]
results1 = solve_ivp(lambda t, x: xdot(t, x, V, kp, ki, R, B0, wg, V_base, XL_base, w_base, 0.5e6/S_base, 0), tspan1, x1, max_step=tstep)
y1 = results1.y

x2 = y1[:,len(y1)]
results2 = solve_ivp(lambda t, x: xdot(t, x, V, kp, ki, R, B0, wg, V_base, XL_base, w_base, 0.5e6/S_base, 200e3/S_base), tspan2, x2, max_step=tstep)
y2 = results2.y

x3 = y2[:,len(y2)]
results3 = solve_ivp(lambda t, x: xdot(t, x, V, kp, ki, R, B0, wg, V_base, XL_base, w_base, 0.25e6/S_base, 200e3/S_base), tspan3, x3, max_step=tstep)
y3 = results3.y

# concatenate results
t0 = results0.t
t1 = results1.t
t2 = results2.t
t3 = results3.t
t = np.concatenate((t0, t1, t2, t3))

theta_g0 = y0[0]
theta_g1 = y1[0]
theta_g2 = y2[0]
theta_g3 = y3[0]
theta_g = np.concatenate((theta_g0, theta_g1, theta_g2, theta_g3))

id0 = y0[1]
id1 = y1[1]
id2 = y2[1]
id3 = y3[1]
i_d = np.concatenate((id0, id1, id2, id3))

iq0 = y0[2]
iq1 = y1[2]
iq2 = y2[2]
iq3 = y3[3]
i_q = np.concatenate((iq0, iq1, iq2, iq3))

# plot results
for i, angle in enumerate(theta_g):
    theta_g[i] = np.mod(angle, 2*pi)

plt.figure(1)
plt.plot(t, theta_g, color='b')
plt.xlabel('Time (sec)')
plt.ylabel('Angle mod 2*pi (rad)')
plt.title('Grid Angle')

plt.figure(2)
plt.plot(t, i_d, label='id',color='b')
plt.plot(t, i_q, label='iq',color='r')
plt.xlabel('Time (sec)')
plt.ylabel('dq-Current (pu)')
plt.legend(('id','iq'))
plt.title("Per Unit Currents")

plt.show()