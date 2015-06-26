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


filename1 = 'strong_scaling/hybrid_data_32_strong_scaling_omp_4096.dat'
#filename2 = 'strong_scaling/hybrid_data_16_strong_scaling_mpi_4096.dat'

data1 = np.loadtxt(filename1) # number of mpi tasks, number of omp threads, time, size, error
#data2 = np.loadtxt(filename2) # number of mpi tasks, number of omp threads, time, size, error

# not actually serial, but with 1 omp thread
serial1 = data1[0,2]
#serial2 = data2[0,2]


if len(data1[0]) > 4: # with error
    plt.errorbar(data1[:,1], serial1/data1[:,2], yerr=data1[:,4])#, label='8 OpenMP threads')
else: # without error
    plt.plot(data1[:,1], serial1/data1[:,2])#, label='8 OpenMP threads')


#if len(data2[0]) > 4: # with error
#    plt.errorbar(data2[:,0], serial2/data2[:,2], yerr=data2[:,4], label='16 OpenMP threads')
#else: # without error
#    plt.plot(data2[:,0], serial2/data2[:,2], label='16 OpenMP threads')
    
    
# annotate with size
#for i in range(0,len(data1)):
#    plt.annotate('%d' % data1[i,3], xy=(data1[i,0],serial1/data2[i,2]))




# plot linear speedup line
plt.plot([0,data1[-1,1]+1], [0,data1[-1,1]+1], label='Linear speedup', color='#BDBDBD', linestyle='--')



plt.ylabel(r'Speedup $t_{1} / t_{n}$')
plt.title('OpenMP strong scaling of hybrid version (N = 4096, transposing with data types)')
plt.legend()
plt.xlim(xmin=0, xmax=data1[-1,1]+1)
plt.ylim(ymax=2)

plt.xlabel('Number of threads')




plt.savefig('strong_scaling/hybrid_strong_scaling_omp_4096.pdf')
plt.show()
