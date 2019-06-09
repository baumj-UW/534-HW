from sim import sim
import matplotlib.pyplot as plt
import numpy as np
import csv

# number of inverters
N = 32
# reference deviation
d = 0.0

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

(t, states, P, Q) = sim(Pref[0], Qref[0])

for n in range(1, N):
    (tn, statesn, Pn, Qn) = sim(Pref[n], Qref[n])
    P = P + Pn
    Q = Q + Qn
    for i, state in enumerate(states):
        if (i > 0):
            states[i] = states[i] + statesn[i]
    
# plot results
plt.figure(1)
plt.plot(t, states[0], color='b')
plt.xlabel('Time (sec)')
plt.ylabel('Angle mod 2*pi (rad)')
plt.title('Grid Angle')
with open ('python_data/grid_angle.csv', mode='w', newline='') as file:
    file_writer = csv.writer(file, delimiter=',')
    for i, time in enumerate(t):
        file_writer.writerow([time, states[0][i]])

plt.figure(2)
plt.plot(t, states[1], label='id',color='b')
plt.plot(t, states[2], label='iq',color='r')
plt.xlabel('Time (sec)')
plt.ylabel('dq-Current (pu)')
plt.legend(('id','iq'))
plt.title("Per Unit Currents")
with open ('python_data/id.csv', mode='w', newline='') as file:
    file_writer = csv.writer(file, delimiter=',')
    for i, time in enumerate(t):
        file_writer.writerow([time, states[1][i]])
with open ('python_data/iq.csv', mode='w', newline='') as file:
    file_writer = csv.writer(file, delimiter=',')
    for i, time in enumerate(t):
        file_writer.writerow([time, states[2][i]])

plt.figure(3)
plt.plot(t, P, label='P', color='b')
plt.plot(t, Q, label='Q', color='r')
plt.xlabel('Time (sec)')
plt.ylabel('Power (pu)')
plt.legend(('P','Q'))
plt.title("Real and Reactive Power")
with open ('python_data/P.csv', mode='w', newline='') as file:
    file_writer = csv.writer(file, delimiter=',')
    for i, time in enumerate(t):
        file_writer.writerow([time, P[i]])
with open ('python_data/Q.csv', mode='w', newline='') as file:
    file_writer = csv.writer(file, delimiter=',')
    for i, time in enumerate(t):
        file_writer.writerow([time, Q[i]])


plt.figure(4)
plt.plot(t, states[5], label='v_c_d',color='b')
plt.plot(t, states[6], label='v_c_q',color='r')
plt.xlabel('Time (sec)')
plt.ylabel('Voltage (pu)')
plt.legend(('v_c_d','v_c_q'))
plt.title("Per Unit LCL Capacitor Voltage")
with open ('python_data/vcd.csv', mode='w', newline='') as file:
    file_writer = csv.writer(file, delimiter=',')
    for i, time in enumerate(t):
        file_writer.writerow([time, states[5][i]])
with open ('python_data/vcq.csv', mode='w', newline='') as file:
    file_writer = csv.writer(file, delimiter=',')
    for i, time in enumerate(t):
        file_writer.writerow([time, states[6][i]])

plt.figure(5)
plt.plot(t, states[7], label='i_g_d',color='b')
plt.plot(t, states[8], label='i_g_q',color='r')
plt.xlabel('Time (sec)')
plt.ylabel('dq-Current (pu)')
plt.legend(('i_g_d','i_g_q'))
plt.title("Per Unit Grid Currents")
with open ('python_data/igd.csv', mode='w', newline='') as file:
    file_writer = csv.writer(file, delimiter=',')
    for i, time in enumerate(t):
        file_writer.writerow([time, states[7][i]])
with open ('python_data/igq.csv', mode='w', newline='') as file:
    file_writer = csv.writer(file, delimiter=',')
    for i, time in enumerate(t):
        file_writer.writerow([time, states[8][i]])

plt.show()