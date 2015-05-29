import numpy as np
import matplotlib.pyplot as plt
import sys

# get filename
if len(sys.argv) < 2:
	filename = raw_input('Data-file name: ')
else:
	filename = sys.argv[1]

#filename2 = 'daint_mpi_datatype.dat'


# load data
data = np.loadtxt(filename) # number of mpi tasks, number of omp threads, time, size
#data2 = np.loadtxt(filename2)

count = len(data[(data[:,0]==2),:])

for i in range(0,len(data)/count):
    plt.plot(data[i*count:(i+1)*count,1], data[i*count:(i+1)*count,2], '-o', label='%d mpi-tasks, N = %d' % (data[i*count,0], data[i*count,3]))

# annotate with size
#for i in range(0,len(data)):
#	plt.annotate('%d' % data[i,3], xy=(data[i,1],data[i,2]))

size = filename.split('_')[-1].split('.')[0]

plt.xlabel('number of OpenMP threads')
plt.ylabel('time')
plt.title('MPI+OpenMP on Piz Daint, N = %s' % size)
plt.legend()
plt.xlim(xmin=0, xmax=data[-1,1]+1) #, xmax=18)

#plt.savefig('daint_hybrid_%s.pdf' % size)
plt.savefig(filename.split('.',1)[0] + '.pdf')
plt.show()
