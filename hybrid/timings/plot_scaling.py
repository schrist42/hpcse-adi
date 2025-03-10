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
data = np.loadtxt(filename) # number of mpi tasks, number of omp threads, time, size, error
#data2 = np.loadtxt(filename2)

count = len(data[(data[:,0]==2),:])

if len(data[0]) > 4:
    for i in range(0,len(data)/count):
        t1 = data[i*count,2]
#        plt.plot(data[i*count:(i+1)*count,1], data[i*count:(i+1)*count,2], '-o', label='%d mpi-tasks, N = %d' % (data[i*count,0], data[i*count,3]))
        plt.errorbar(data[i*count:(i+1)*count,1], t1/data[i*count:(i+1)*count,2], yerr=data[i*count:(i+1)*count,4], label='%d MPI processes, N = %d' % (data[i*count,0], data[i*count,3]))
else:
    for i in range(0,len(data)/count):
        t1 = data[i*count,2]
        plt.plot(data[i*count:(i+1)*count-1,1], t1/data[i*count:(i+1)*count-1,2], '-o', label='%d mpi-tasks, N = %d' % (data[i*count,0], data[i*count,3]))
#        plt.errorbar(data[i*count:(i+1)*count,1], data[i*count:(i+1)*count,2], yerr=data[i*count:(i+1)*count,4], label='%d mpi-tasks, N = %d' % (data[i*count,0], data[i*count,3]))

# annotate with size
#for i in range(0,len(data)):
#    plt.annotate('%d' % data[i,3], xy=(data[i,1],data[i,2]))

print data[count,1]
plt.plot([0,data[-1,1]+1], [0,data[-1,1]+1], label='Linear speedup', color='#BDBDBD', linestyle='--')

size = filename.split('_')[-1].split('.')[0]

plt.xlabel('Number of threads')
plt.ylabel(r'Speedup $t_1 / t_n$')
plt.title('OpenMP strong scaling of hybrid version\n(N = %s, transposing blocks locally)' % size)
plt.legend(loc='upper left')
plt.xlim(xmin=0, xmax=data[-1,1]) #, xmax=18)
plt.ylim(ymax=3.5)

#plt.savefig('daint_hybrid_%s.pdf' % size)
plt.savefig(filename.split('.',1)[0] + '.pdf')
plt.show()
