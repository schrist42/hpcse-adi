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


filename1 = 'strong_scaling/hybrid_data_8_strong_scaling_mpi_512.dat'
filename2 = 'strong_scaling/hybrid_data_8_strong_scaling_mpi_1024.dat'
filename3 = 'strong_scaling/hybrid_data_8_strong_scaling_mpi_2048.dat'
filename4 = 'strong_scaling/hybrid_data_8_strong_scaling_mpi_4096.dat'

data1 = np.loadtxt(filename1) # number of mpi tasks, number of omp threads, time, size, error
data2 = np.loadtxt(filename2) 
data3 = np.loadtxt(filename3)
data4 = np.loadtxt(filename4)

# not actually serial, but with 1 mpi task
serial1 = data1[0,2]
serial2 = data2[0,2]
serial3 = data3[0,2]
serial4 = data4[0,2]


if len(data1[0]) > 4: # with error
    plt.errorbar(data1[:,0], serial1/data1[:,2], yerr=data1[:,4], label='N = 512')
    plt.errorbar(data2[:,0], serial2/data2[:,2], yerr=data2[:,4], label='N = 1024')
    plt.errorbar(data3[:,0], serial3/data3[:,2], yerr=data3[:,4], label='N = 2048')
    plt.errorbar(data4[:,0], serial4/data4[:,2], yerr=data4[:,4], label='N = 4096')
else: # without error
    plt.plot(data1[:,0], serial1/data1[:,2], label='8 OpenMP threads')


#if len(data2[0]) > 4: # with error
#    plt.errorbar(data2[:,0], serial2/data2[:,2], yerr=data2[:,4], label='16 OpenMP threads')
#else: # without error
#    plt.plot(data2[:,0], serial2/data2[:,2], label='16 OpenMP threads')
    
    
# annotate with size
#for i in range(0,len(data1)):
#    plt.annotate('%d' % data1[i,3], xy=(data1[i,0],serial1/data2[i,2]))




# plot linear speedup line
plt.plot([0,data1[-1,0]+1], [0,data1[-1,0]+1], label='Linear speedup', color='#BDBDBD', linestyle='--')



plt.ylabel(r'Speedup $t_{1} / t_{n}$')
plt.title('MPI strong scaling of hybrid version\n(8 OpenMP threads, transposing with data types)')
plt.legend()
plt.xlim(xmin=0, xmax=data1[-1,0])

plt.xlabel('Number of processes')




plt.savefig('strong_scaling/hybrid_strong_scaling_mpi_sizes.pdf')
plt.show()
