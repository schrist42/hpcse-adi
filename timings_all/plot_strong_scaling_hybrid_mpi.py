import numpy as np
import matplotlib.pyplot as plt
import sys

# get filename
#if len(sys.argv) < 2:
#    filename = raw_input('Data-file name: ')
#else:
#    filename = sys.argv[1]


#if filename.startswith('strong_scaling/'):
#    filename = filename.split('/',1)[1]
#filename = 'strong_scaling/' + filename
#filename = filename.split('.dat',1)[0] # remove ending

#nthreads = 16
size = 4096

# load data
data_serial = np.loadtxt('strong_scaling/omp_strong_scaling_4096.dat') # threads, time, size, error

ind_serial = np.where(data_serial[:,0] == 8)[0]
t1 = data_serial[ind_serial,1][0]

ind_serial = np.where(data_serial[:,0] == 16)[0]
t2 = data_serial[ind_serial,1][0]


data1 = np.loadtxt('strong_scaling/hybrid_lt_8_strong_scaling_mpi_4096.dat') # number of mpi tasks, number of omp threads, time, size, error
data2 = np.loadtxt('strong_scaling/hybrid_lt_16_strong_scaling_mpi_4096.dat') # number of mpi tasks, number of omp threads, time, size, error


if len(data1[0]) > 4: # with error
    plt.errorbar(data1[:,0], t1/data1[:,2], yerr=data1[:,4], label='8 OpenMP threads')
else: # without error
    plt.plot(data1[:,0], t1/data1[:,2], label='8 OpenMP threads')


if len(data2[0]) > 4: # with error
    plt.errorbar(data2[:,0], t2/data2[:,2], yerr=data2[:,4], label='16 OpenMP threads')
else: # without error
    plt.plot(data2[:,0], t2/data2[:,2], label='16 OpenMP threads')
    
    
# annotate with size
for i in range(0,len(data1)):
    plt.annotate('%d' % data1[i,3], xy=(data1[i,0],t1/data1[i,2]))




# plot linear speedup line
plt.plot([0,data1[-1,0]+1], [0,data1[-1,0]+1], label='Linear speedup', color='#BDBDBD', linestyle='--')



plt.ylabel(r'Speedup $t_{1} / t_{parallel}$ per step')
plt.title('Strong scaling/speedup of hybrid version (different numbers of MPI tasks) on Piz Daint (N = 4096)')
plt.legend()
plt.xlim(xmin=0, xmax=data1[-1,0]+1)

plt.xlabel('Number of tasks')




plt.savefig('hybrid_strong_scaling_mpi_4096.pdf')
plt.show()
