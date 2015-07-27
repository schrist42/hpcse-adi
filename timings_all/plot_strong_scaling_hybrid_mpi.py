import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
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


filename1 = 'strong_scaling/hybrid_data_8_strong_scaling_mpi_4096.dat'
filename2 = 'strong_scaling/hybrid_data_16_strong_scaling_mpi_4096.dat'

data1 = np.loadtxt(filename1) # number of mpi tasks, number of omp threads, time, size, error
data2 = np.loadtxt(filename2) # number of mpi tasks, number of omp threads, time, size, error

# not actually serial, but with 1 mpi task
serial1 = data1[0,2]
serial2 = data2[0,2]

fig1, ax1 = plt.subplots()

if len(data1[0]) > 4: # with error
    ax1.errorbar(data1[:,0], serial1/data1[:,2], yerr=data1[:,4], label='8 OpenMP threads')
else: # without error
    ax1.plot(data1[:,0], serial1/data1[:,2], label='8 OpenMP threads')


if len(data2[0]) > 4: # with error
    ax1.errorbar(data2[:,0], serial2/data2[:,2], yerr=data2[:,4], label='16 OpenMP threads')
else: # without error
    ax1.plot(data2[:,0], serial2/data2[:,2], label='16 OpenMP threads')
    
    
# annotate with size
#for i in range(0,len(data1)):
#    ax1.annotate('%d' % data1[i,3], xy=(data1[i,0],serial1/data2[i,2]))




# plot linear speedup line
ax1.plot([0,data1[-1,0]+1], [0,data1[-1,0]+1], label='Linear speedup', color='#BDBDBD', linestyle='--')



ax1.set_ylabel(r'Speedup $t_{1} / t_{n}$')
ax1.set_title('MPI strong scaling of hybrid version\n(N = 4096, transposing with data types)')
plt.legend()
ax1.set_xlim(xmin=0, xmax=data1[-1,0]+1)

ax1.set_xlabel('Number of processes')

labels = [''] * len(data1[:,0])
labels[0] = '%d' % int(data1[0,0])
labels[-4] = '%d' % int(data1[-4,0])
labels[-3] = '%d' % int(data1[-3,0])
labels[-2] = '%d' % int(data1[-2,0])
labels[-1] = '%d' % int(data1[-1,0])

ax1.xaxis.set_major_locator(ticker.FixedLocator(data1[:,0]))
ax1.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
ax1.xaxis.set_minor_locator(ticker.NullLocator())

plt.savefig('strong_scaling/hybrid_strong_scaling_mpi_4096.pdf')
plt.show()
