import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import sys

# get filename
#if len(sys.argv) < 2:
#    filename = raw_input('Data-file name: ')
#else:
#    filename = sys.argv[1]
#
#
#if filename.startswith('serial/'):
#    filename = filename.split('/',1)[1]
#filename = 'serial/' + filename
#filename = filename.split('.dat',1)[0] # remove ending
#
#parallel_type = filename.split('/')[1].split('_')[0]

# load data
data_lt = np.loadtxt('mpi/mpi_lt_32_size.dat') # tasks, time, size, error
data_datatype = np.loadtxt('mpi/mpi_data_32_size.dat') # tasks, time, size, error

fig1, ax1 = plt.subplots()

if len(data_lt[0]) > 3: # with error
    ax1.errorbar(data_lt[:,2], data_lt[:,1], yerr=data_lt[:,3], label='block-wise local transpose')
    ax1.errorbar(data_datatype[:,2], data_datatype[:,1], yerr=data_datatype[:,3], label='different send/receive data types')
else: # without error
    plt.loglog(data_lt[:,2], data_lt[:,1], 'x-', label='block-wise local transpose')
    plt.loglog(data_datatype[:,2], data_datatype[:,1], 'x-', label='different send/receive data types')

# annotate with size
#for i in range(0,len(data)):
#    plt.annotate('%d' % data[i,2], xy=(data[i,0],t1/data[i,1]))

#size = filename.split('_')[-1].split('.',1)[0]

logx = np.log(data_lt[:,2])
logy = np.log(data_lt[:,1])
coeffs = np.polyfit(logx,logy,deg=1)
print "Scaling for local transpose: %f" % coeffs[0]

logx = np.log(data_datatype[:,2])
logy = np.log(data_datatype[:,1])
coeffs = np.polyfit(logx,logy,deg=1)
print "Scaling for data types: %f" % coeffs[0]

ax1.set_xlabel('Size')
ax1.set_ylabel('Time per step')
ax1.set_title('System size scaling of MPI version\n(32 processes)')
ax1.legend(loc='upper left')
ax1.set_xscale('log')
ax1.set_yscale('log')
#ax1.set_xticks([], minor=False)
#ax1.set_xticks(data_lt[:,2], minor=False) #[data_lt[0,2], data_lt[-1,2]])
#
#labels = [''] * len(ax1.get_xticklabels())
#labels[0] = int(data_lt[0,2])
#labels[-1] = int(data_lt[-1,2])
#
#ax1.set_xticklabels(labels)


labels = [''] * len(data_lt[:,2])
labels[0] = '%d' % int(data_lt[0,2])
labels[1] = '%d' % int(data_lt[1,2])
labels[2] = '%d' % int(data_lt[2,2])
labels[4] = '%d' % int(data_lt[4,2])
labels[-1] = '%d' % int(data_lt[-1,2])

ax1.xaxis.set_major_locator(ticker.FixedLocator(data_lt[:,2]))
ax1.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
ax1.xaxis.set_minor_locator(ticker.NullLocator())

ax1.set_xlim(xmin=0, xmax=max(data_lt[:,2]))

plt.savefig('mpi/mpi_transpose_32.pdf')
plt.show()
