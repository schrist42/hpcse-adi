import numpy as np
import matplotlib.pyplot as plt
import sys

# get filename
if len(sys.argv) < 2:
        filename = raw_input('Data-file name: ')
else:
        filename = sys.argv[1]


# load data
data = np.loadtxt(filename) # number of mpi tasks, time, size


plt.plot(data[:,0], data[:,1], '-o')#, label='%d mpi-tasks, N = %d' % (data[i*count,0], data[i*count,3]))

# annotate with size
for i in range(0,len(data)):
        plt.annotate('%d' % data[i,2], xy=(data[i,0],data[i,1]))

size = filename.split('_')[-1].split('.',1)[0]

plt.xlabel('number of CPUs')
plt.ylabel('time')
plt.title('MPI on Piz Daint, N = %s' % size)
#plt.legend()
plt.xlim(xmin=0, xmax=data[-1,0]+1) #, xmax=18)

#plt.savefig('daint_hybrid_%s.pdf' % size)
plt.savefig(filename.split('.',1)[0] + '.pdf')
plt.show()
