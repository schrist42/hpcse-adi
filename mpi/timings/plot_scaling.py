import numpy as np
import matplotlib.pyplot as plt
import sys

# get filename
#if len(sys.argv) < 2:
#	filename = raw_input('Data-file name: ')
#else:
#	filename = sys.argv[1]

filename1 = 'daint_mpi_localtranspose_256.dat'
filename2 = 'daint_mpi_localtranspose_512.dat'
#filename2 = 'daint_mpi_datatype.dat'


# load data
data1 = np.loadtxt(filename1) # number of cpus, time, size
data2 = np.loadtxt(filename2) # number of cpus, time, size
#data2 = np.loadtxt(filename2)

plt.plot(data1[:,0], data1[:,1], '-o', label='N = %d' % data1[0,2])
plt.plot(data2[:,0], data2[:,1], '-o', label='N = %d' % data2[0,2])
#plt.plot(data2[:,0], data2[:,1], 'o', label='Datatype')

# annotate with size
for i in range(0,len(data1)):
	plt.annotate('%d' % data1[i,2], xy=(data1[i,0],data1[i,1]))
	plt.annotate('%d' % data2[i,2], xy=(data2[i,0],data2[i,1]))


plt.xlabel('number of CPUs')
plt.ylabel('time')
plt.title('MPI on Piz Daint')
plt.legend()
plt.xlim(xmin=0, xmax=18)

plt.savefig('daint_mpi_sizes.pdf')
plt.show()
