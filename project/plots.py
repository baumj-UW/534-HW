import csv
import os
import matplotlib.pyplot as plt

t_plecs = list()
data_plecs = list()

for filename in os.listdir('plecs_data'):
    print(filename)
    data_plecs.append(list())

for i, filename in enumerate(os.listdir('plecs_data')):
    with open('plecs_data/' + filename) as file:
        file_reader = csv.reader(file, delimiter=',')
        line = 0 # skip first line
        for row in file_reader:
            if line == 0:
                line += 1
            else:
                if i == 0:
                    t_plecs.append(float(row[0]))
                data_plecs[i].append(float(row[1]))

t_py = list()
data_py = list()

for filename in os.listdir('python_data'):
    print(filename)
    data_py.append(list())

for i, filename in enumerate(os.listdir('python_data')):
    with open('python_data/' + filename) as file:
        file_reader = csv.reader(file, delimiter=',')
        line = 0 # skip first line
        for row in file_reader:
            if line == 0:
                line += 1
            else:
                if i == 0:
                    t_py.append(float(row[0]))
                data_py[i].append(float(row[1]))


# plot results
plt.figure(1)
plt.plot(t_py, data_py[0], label='id',color='b')
plt.plot(t_py, data_py[3], label='iq',color='b')
plt.plot(t_plecs, data_plecs[0], label='id',color='r')
plt.plot(t_plecs, data_plecs[3], label='iq',color='r')
plt.xlabel('Time (sec)')
plt.ylabel('dq-Current (pu)')
# plt.legend(('id','iq'))
plt.title("Per Unit Currents")

plt.figure(2)
plt.plot(t_py, data_py[4], label='P', color='b')
plt.plot(t_py, data_py[5], label='Q', color='b')
plt.plot(t_plecs, data_plecs[4], label='P', color='r')
plt.plot(t_plecs, data_plecs[5], label='Q', color='r')
plt.xlabel('Time (sec)')
plt.ylabel('Power (pu)')
# plt.legend(('P','Q'))
plt.title("Real and Reactive Power")

plt.figure(3)
plt.plot(t_py, data_py[1], label='i_g_d',color='b')
plt.plot(t_py, data_py[2], label='i_g_q',color='b')
plt.plot(t_plecs, data_plecs[1], label='i_g_d',color='r')
plt.plot(t_plecs, data_plecs[2], label='i_g_q',color='r')
plt.xlabel('Time (sec)')
plt.ylabel('dq-Current (pu)')
# plt.legend(('i_g_d','i_g_q'))
plt.title("Per Unit Grid Currents")

plt.show()