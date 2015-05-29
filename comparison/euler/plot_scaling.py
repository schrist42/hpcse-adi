import numpy as np
import matplotlib.pyplot as plt
import sys

# get filename
#if len(sys.argv) < 2:
#	filename = raw_input('Data-file name: ')
#else:
#	filename = sys.argv[1]

filename1 = 'euler_omp_256.dat'
filename2 = 'euler_mpi_localtranspose_256.dat'
#filename3 = 'euler_hybrid_256.dat'


# load data
data1 = np.loadtxt(filename1) # number of threads, time, size
data2 = np.loadtxt(filename2) # number of cpus, time, size
#data3 = np.loadtxt(filename3)

plt.plot(data1[:,0], data1[:,1], '-o', label='OMP')
plt.plot(data2[:,0], data2[:,1], '-o', label='MPI')
#plt.plot(data3[:,0], data3[:,1], 'o', label='Hybrid')

# annotate with size
for i in range(0,len(data2)):
	#plt.annotate('%d' % data1[i,2], xy=(data1[i,0],data1[i,1]))
	plt.annotate('%d' % data2[i,2], xy=(data2[i,0],data2[i,1]))


plt.xlabel('number of OMP threads or MPI processes')
plt.ylabel('time')
plt.title('Comparison on Euler')
plt.legend()
plt.xlim(xmin=0, xmax=25)

plt.savefig('comparison_256.pdf')
plt.show()
