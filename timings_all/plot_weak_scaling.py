import numpy as np
import matplotlib.pyplot as plt
import sys

# get filename
if len(sys.argv) < 2:
        filename = raw_input('Data-file name: ')
else:
        filename = sys.argv[1]


if filename.startswith('weak_scaling/'):
    filename = filename.split('/',1)[1]
filename = 'weak_scaling/' + filename
filename = filename.split('.dat',1)[0] # remove ending

parallel_type = filename.split('/')[1].split('_')[0]

# load data
data = np.loadtxt(filename + '.dat') # number of mpi tasks, time, size, error

if len(data[0]) > 3: # with error
    plt.errorbar(data[:,0], data[:,1], yerr=data[:,3])#, label='%d mpi-tasks, N = %d' % (data[i*count,0], data[i*count,3]))
else: # without error
    plt.plot(data[:,0], data[:,1])#, label='%d mpi-tasks, N = %d' % (data[i*count,0], data[i*count,3]))

# annotate with size
for i in range(0,len(data)):
    plt.annotate('%d' % data[i,2], xy=(data[i,0],data[i,1]))

size = filename.split('_')[-1].split('.',1)[0]

plt.xlabel('x times threads or original size (N=1024)')
plt.ylabel('time per step')
plt.title('Weak scaling of %s on Piz Daint' % parallel_type)
#plt.legend()
#plt.xlim(xmin=0, xmax=data[-1,0]+1) #, xmax=18)

plt.savefig(filename + '.pdf')
plt.show()
