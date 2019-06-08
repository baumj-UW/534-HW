from constants import *
from xdot import xdot

import numpy as np
from scipy.integrate import solve_ivp


def sim(Pref, Qref):
    """
    Simulates an inverter from t=0 to t=20*Tg in four intervals,
    each of which are 5*Tg in length
    
    Arguments:
        Pref {list} -- list containing 4 reference real powers, for each interval
        Qref {list} -- list containing 4 reference imaginary powers, for each interval
    
    Returns:
        list -- list containing the arrays (t, states, P, Q)
    """

    #unpack Pref and Qref
    (Pref0, Pref1, Pref2, Pref3) = Pref
    (Qref0, Qref1, Qref2, Qref3) = Qref

    # maximum step size for ode solver
    tstep = 5*Tg/50 #1e3

    # time intervals
    tspan0 = [0, 5*Tg]
    tspan1 = [5*Tg, 10*Tg]
    tspan2 = [10*Tg, 15*Tg]
    tspan3 = [15*Tg, 20*Tg]

    # initial values
    # thetag, id, iq, zd, zq, vcd, vcq, igd, igq
    x0 = 3.14159265e+01, 1.41421356e+00, 0.0, 2.14832789e-02, 0.0, 3.89241188e-01, -3.54671089e-01, 1.26684478e+00, -1.62349027e-01

    # simulate!
    results0 = solve_ivp(lambda t, x: xdot(t, x, V, kp, ki, R, Rg, B0, R0, wg, V_base, XL, XLg, w_base, Pref0, Qref0), tspan0, x0, max_step=tstep)
    y0 = results0.y

    x1 = y0[:,y0.shape[1]-1]
    results1 = solve_ivp(lambda t, x: xdot(t, x, V, kp, ki, R, Rg, B0, R0, wg, V_base, XL, XLg, w_base, Pref1, Qref1), tspan1, x1, max_step=tstep)
    y1 = results1.y

    x2 = y1[:,y1.shape[1]-1]
    results2 = solve_ivp(lambda t, x: xdot(t, x, V, kp, ki, R, Rg, B0, R0, wg, V_base, XL, XLg, w_base, Pref2, Qref2), tspan2, x2, max_step=tstep)
    y2 = results2.y

    x3 = y2[:,y2.shape[1]-1]
    results3 = solve_ivp(lambda t, x: xdot(t, x, V, kp, ki, R, Rg, B0, R0, wg, V_base, XL, XLg, w_base, Pref3, Qref2), tspan3, x3, max_step=tstep)
    y3 = results3.y

    # concatenate results
    t = np.concatenate((results0.t, results1.t, results2.t, results3.t))
    states = list()
    for i in range(0,len(y0)):
        data = np.concatenate((y0[i], y1[i], y2[i], y3[i]))
        states.append(data)
    
    for i, angle in enumerate(states[0]):
        states[0][i] = np.mod(angle, 2*pi)
    
    P = list()
    Q = list()

    for i, cur in enumerate(states[0]):
        P.append((V / V_base) * states[7][i])
        Q.append((V / V_base) * states[8][i])
    
    return (t, states, P, Q)
