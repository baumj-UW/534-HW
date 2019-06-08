from sim import sim
import matplotlib.pyplot as plt
import numpy as np

# number of inverters
N = 10
# reference deviation
d = 0.1

Pref = list()
Qref = list()
for n in range(0, N):
    # random factor between (1-d) and (1+d)
    factor = ((1.0 - d) + np.random.rand() * (2 * d))
    Pref0 = 1.0 * factor
    Pref1 = 0.5 * factor
    Pref2 = 0.5 * factor
    Pref3 = 0.25 * factor

    Qref0 = 0.0 * factor
    Qref1 = 0.0 * factor
    Qref2 = 0.2 * factor
    Qref3 = 0.2 * factor
    Pref.append((Pref0, Pref1, Pref2, Pref3))
    Qref.append((Qref0, Qref1, Qref2, Qref3))

t = list()
states = list()
P = list()
Q = list()
for n in range(0, N):
    (tn, statesn, Pn, Qn) = sim(Pref[n], Qref[n])
    t.append(tn)
    states.append(statesn)
    P.append(Pn)
    Q.append(Qn)
    
# plot results
plt.figure(1)
for n in range(0, N):
    plt.plot(t[n], states[n][0], color='b')
plt.xlabel('Time (sec)')
plt.ylabel('Angle mod 2*pi (rad)')
plt.title('Grid Angle')

plt.figure(2)
for n in range(0, N):
    plt.plot(t[n], states[n][1], label='id',color='b')
    plt.plot(t[n], states[n][2], label='iq',color='r')
plt.xlabel('Time (sec)')
plt.ylabel('dq-Current (pu)')
plt.legend(('id','iq'))
plt.title("Per Unit Currents")

plt.figure(3)
for n in range(0, N):
    plt.plot(t[n], P[n], label='P', color='b')
    plt.plot(t[n], Q[n], label='Q', color='r')
plt.xlabel('Time (sec)')
plt.ylabel('Power (pu)')
plt.legend(('P','Q'))
plt.title("Real and Reactive Power")

plt.figure(4)
for n in range(0, N):
    plt.plot(t[n], states[n][5], label='v_c_d',color='b')
    plt.plot(t[n], states[n][6], label='v_c_q',color='r')
plt.xlabel('Time (sec)')
plt.ylabel('Voltage (pu)')
plt.legend(('v_c_d','v_c_q'))
plt.title("Per Unit LCL Capacitor Voltage")

plt.figure(5)
for n in range(0, N):
    plt.plot(t[n], states[n][7], label='i_g_d',color='b')
    plt.plot(t[n], states[n][8], label='i_g_q',color='r')
plt.xlabel('Time (sec)')
plt.ylabel('dq-Current (pu)')
plt.legend(('i_g_d','i_g_q'))
plt.title("Per Unit Grid Currents")

plt.show()