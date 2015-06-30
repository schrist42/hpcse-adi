import numpy as np
import matplotlib.pyplot as plt
import sys

# get filename
filename = 'weak_scaling/omp_weak_scaling.dat'

# load data
data_serial = np.loadtxt('serial/serial_vec.dat') # time, size, error
t1 = data_serial[1,0]
data = np.loadtxt(filename) # number of omp threads, time, size, error
#t1 = data[0,1]

if len(data[0]) > 3: # with error
    plt.errorbar(data[:,0], t1/data[:,1], yerr=data[:,3])#, label='%d mpi-tasks, N = %d' % (data[i*count,0], data[i*count,3]))
else: # without error
    plt.plot(data[:,0], t1/data[:,1])#, label='%d mpi-tasks, N = %d' % (data[i*count,0], data[i*count,3]))

# annotate with size
for i in range(0,len(data)):
    plt.annotate('%d' % data[i,2], xy=(data[i,0],t1/data[i,1]))

size = filename.split('_')[-1].split('.',1)[0]

# plot linear speedup line
plt.plot([0,data[-1,0]], [1,1], label='Linear weak scaling', color='#BDBDBD', linestyle='--')

plt.xlabel('x times threads or original size (N=1024)')
plt.ylabel(r'$t_{1} / t_{n}$ per step')
plt.title('Weak scaling of OpenMP version') # % parallel_type)
#plt.legend()
plt.xlim(xmin=0, xmax=data[-1,0]) #, xmax=18)

plt.savefig('weak_scaling/omp_weak_scaling.pdf')
plt.show()
