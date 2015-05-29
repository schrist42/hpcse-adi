import numpy as np
import matplotlib.pyplot as plt
import sys

# get filename
#if len(sys.argv) < 2:
#	filename = raw_input('Data-file name: ')
#else:
#	filename = sys.argv[1]

filename1 = 'euler_omp_256.dat'
filename2 = 'euler_omp_384.dat'
filename3 = 'euler_omp_512.dat'

# load data
data1 = np.loadtxt(filename1) # number of cpus, time, size
data2 = np.loadtxt(filename2)
data3 = np.loadtxt(filename3)
plt.plot(data1[:,0], data1[:,1], 'o', label='N256')
plt.plot(data2[:,0], data2[:,1], 'o', label='N384')
plt.plot(data3[:,0], data3[:,1], 'o', label='N512')

# annotate with size
#for i in range(0,len(data1)):
#	plt.annotate('%d' % data1[i,2], xy=(data1[i,0],data1[i,1]))


plt.xlabel('Number of threads')
plt.ylabel('Time [s]')
plt.title('OMP on EULER')
plt.legend()

plt.savefig('euler_omp_sizes.pdf')
plt.show()
