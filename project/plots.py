import csv
import os
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

t_plecs = list()
data_plecs = list()

for filename in os.listdir('plecs_data'):
    # print(filename)
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

t_plecs_scaled = list()
data_plecs_scaled = list()

for filename in os.listdir('plecs_data_scaled'):
    # print(filename)
    data_plecs_scaled.append(list())

for i, filename in enumerate(os.listdir('plecs_data_scaled')):
    with open('plecs_data_scaled/' + filename) as file:
        file_reader = csv.reader(file, delimiter=',')
        line = 0 # skip first line
        for row in file_reader:
            if line == 0:
                line += 1
            else:
                if i == 0:
                    t_plecs_scaled.append(float(row[0]))
                data_plecs_scaled[i].append(float(row[1]))

t_py = list()
data_py = list()

for filename in os.listdir('python_data'):
    # print(filename)
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

t_py_scaled = list()
data_py_scaled = list()

for filename in os.listdir('python_data_scaled'):
    # print(filename)
    data_py_scaled.append(list())

for i, filename in enumerate(os.listdir('python_data_scaled')):
    with open('python_data_scaled/' + filename) as file:
        file_reader = csv.reader(file, delimiter=',')
        line = 0 # skip first line
        for row in file_reader:
            if line == 0:
                line += 1
            else:
                if i == 0:
                    t_py_scaled.append(float(row[0]))
                data_py_scaled[i].append(float(row[1]))


# plot results
plt.figure(1)
plt.plot(t_py, data_py[0], label='i, Python Aggregate Model', color='C0')
plt.plot(t_py, data_py[3], color='C0')
plt.plot(t_py_scaled, data_py_scaled[0], label='i, Python Scaled Model',color='C1')
plt.plot(t_py_scaled, data_py_scaled[3], color='C1')
plt.plot(t_plecs, data_plecs[0], label='i, PLECS Aggregated Circuit Equivalent',color='C2', linewidth='2.0')
plt.plot(t_plecs, data_plecs[3], color='C2', linewidth='2.0')
plt.plot(t_plecs_scaled, data_plecs_scaled[0], label='i, PLECS Scaled Circuit Equivalent',color='C3')
plt.plot(t_plecs_scaled, data_plecs_scaled[3], color='C3')
plt.xlabel('Time (sec)')
plt.ylabel('dq-Current (pu)')
plt.legend(('Aggregate Model (Python)','Aggregate Circuit Equivalent (PLECS)', 'Scaled Model (Python)', 'Scaled Circuit Equivalent (PLECS)'))
plt.xlim((0, 20/60))
plt.title("Per Unit Currents")

plt.figure(2)
plt.plot(t_py, data_py[4], label='P/Q, Python Aggregate Model', color='C0')
plt.plot(t_py, data_py[5], color='C0')
plt.plot(t_py_scaled, data_py_scaled[4], label='P/Q, Python Scaled Model',color='C1')
plt.plot(t_py_scaled, data_py_scaled[5], color='C1')
plt.plot(t_plecs, data_plecs[4], label='P/Q, PLECS Aggregated Circuit Equivalent',color='C2', linewidth='2.0')
plt.plot(t_plecs, data_plecs[5], color='C2', linewidth='2.0')
plt.plot(t_plecs_scaled, data_plecs_scaled[4], label='P/Q, PLECS Scaled Circuit Equivalent',color='C3')
plt.plot(t_plecs_scaled, data_plecs_scaled[5], color='C3')
plt.xlabel('Time (sec)')
plt.ylabel('Power (pu)')
plt.legend()
plt.xlim((0, 20/60))
plt.title("Real and Reactive Power")

plt.figure(3)
plt.plot(t_py, data_py[1], label='ig, Python Aggregate Model', color='C0')
plt.plot(t_py, data_py[2], color='C0')
plt.plot(t_py_scaled, data_py_scaled[1], label='ig, Python Scaled Model',color='C1')
plt.plot(t_py_scaled, data_py_scaled[2], color='C1')
plt.plot(t_plecs, data_plecs[1], label='ig, PLECS Aggregated Circuit Equivalent',color='C2', linewidth='2.0')
plt.plot(t_plecs, data_plecs[2], color='C2', linewidth='2.0')
plt.plot(t_plecs_scaled, data_plecs_scaled[1], label='ig, PLECS Scaled Circuit Equivalent',color='C3')
plt.plot(t_plecs_scaled, data_plecs_scaled[2], color='C3')
plt.xlabel('Time (sec)')
plt.ylabel('dq-Current (pu)')
plt.legend()
plt.xlim((0, 20/60))
plt.title("Per Unit Grid Currents")

plt.show()