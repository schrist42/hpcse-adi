import numpy as np
import matplotlib.pyplot as plt
import sys

# get filename
#if len(sys.argv) < 2:
#	filename = raw_input('Data-file name: ')
#else:
#	filename = sys.argv[1]

filename1 = 'daint_hybrid_localtranspose_512.dat'
#filename2 = 'daint_mpi_datatype.dat'


# load data
data = np.loadtxt(filename1) # number of mpi tasks, number of omp threads, time, size
#data2 = np.loadtxt(filename2)

count = len(data[(data[:,0]==2),:])

for i in range(0,len(data)/count):
    plt.plot(data[i*count:(i+1)*count,1], data[i*count:(i+1)*count,2], '-o', label='%d mpi-tasks' % data[i*count,0])

# annotate with size
for i in range(0,len(data)):
	plt.annotate('%d' % data[i,3], xy=(data[i,1],data[i,2]))


plt.xlabel('number of OpenMP threads')
plt.ylabel('time')
plt.title('MPI+OpenMP on Piz Daint, N = 512')
plt.legend()
#plt.xlim(xmin=0, xmax=18)

plt.savefig('daint_hybrid_512.pdf')
plt.show()
