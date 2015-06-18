import numpy as np
import matplotlib.pyplot as plt
import sys

# get filename
if len(sys.argv) < 2:
    filename = raw_input('Data-file name: ')
else:
    filename = sys.argv[1]


if filename.startswith('strong_scaling/'):
    filename = filename.split('/',1)[1]
filename = 'strong_scaling/' + filename
filename = filename.split('.dat',1)[0] # remove ending

nthreads = 16
#parallel_type = filename.split('/')[1].split('_')[0]
size = filename.split('_')[-1]
ptype_hybrid = filename.split('_')[-2]

# load data
data_serial = np.loadtxt('strong_scaling/omp_strong_scaling_%s.dat' % size) # time, size, error
ind_serial = np.where(data_serial[:,0] == nthreads)[0]
#t1 = data_serial.where(data_serial[:,1]=4096)
t1 = data_serial[ind_serial,1][0]

data = np.loadtxt(filename + '.dat') # number of mpi tasks, number of omp threads, time, size, error

if ptype_hybrid == 'omp':
    if len(data[0]) > 4: # with error
        plt.errorbar(data[:,1], t1/data[:,2], yerr=data[:,4])#, label='%d mpi-tasks, N = %d' % (data[i*count,0], data[i*count,3]))
    else: # without error
        plt.plot(data[:,1], t1/data[:,2])#, label='%d mpi-tasks, N = %d' % (data[i*count,0], data[i*count,3]))
else:
    if len(data[0]) > 4: # with error
        plt.errorbar(data[:,0], t1/data[:,2], yerr=data[:,4])#, label='%d mpi-tasks, N = %d' % (data[i*count,0], data[i*count,3]))
    else: # without error
        plt.plot(data[:,0], t1/data[:,2])#, label='%d mpi-tasks, N = %d' % (data[i*count,0], data[i*count,3]))
    # annotate with size
    for i in range(0,len(data)):
        plt.annotate('%d' % data[i,3], xy=(data[i,0],t1/data[i,2]))


# plot linear speedup line
plt.plot([0,data[-1,0]+1], [0,data[-1,0]+1], label='Linear speedup', color='#BDBDBD', linestyle='--')


#size = filename.split('_')[-1].split('.',1)[0]

plt.ylabel(r'Speedup $t_{serial} / t_{parallel}$ per step')
plt.title('Strong scaling/speedup of hybrid version (different numbers of MPI tasks for %d OpenMP threads) on Piz Daint' % nthreads)
#plt.legend()
plt.xlim(xmin=0, xmax=data[-1,0]+1) #, xmax=18)
#plt.ylim(ymax=data[-1,1]+1) #, xmax=18)


if ptype_hybrid == 'omp':
    plt.xlabel('Number of threads')
elif ptype_hybrid == 'mpi':
    plt.xlabel('Number of tasks')




plt.savefig(filename + '.pdf')
plt.show()
