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
data1 = np.loadtxt('omp/omp_size_scaling_8.dat') # threads, time, size, error
data2 = np.loadtxt('mpi/mpi_data_32_size.dat') # tasks, time, size, error
data3 = np.loadtxt('hybrid/hybrid_size_data.dat') # ntasks, nthreads, time, size, err

fig1, ax1 = plt.subplots()

if len(data1[0]) > 2: # with error
    ax1.errorbar(data1[:,2], data1[:,1], yerr=data1[:,3], label='OpenMP version, 8 threads')
    ax1.errorbar(data2[:,2], data2[:,1], yerr=data2[:,3], label='MPI version, 32 processes')
    ax1.errorbar(data3[:,3], data3[:,2], yerr=data3[:,4], label='Hybrid version, 32 processes, 8 threads, transpose with data type')
else: # without error
    ax1.plot(data1[:,1], data1[:,0], label='vectorized (-ftree-vectorize)')

# annotate with size
#for i in range(0,len(data)):
#    plt.annotate('%d' % data[i,2], xy=(data[i,0],t1/data[i,1]))

#size = filename.split('_')[-1].split('.',1)[0]

logx = np.log(data1[:,2])
logy = np.log(data1[:,1])
coeffs = np.polyfit(logx,logy,deg=1)
poly = np.poly1d(coeffs)

yfit = lambda x: np.exp(poly(np.log(x)))

print "Scaling for OpenMP version: %f" % coeffs[0]
#fitlabel = r'$\alpha =$%f' % coeffs[0]
#ax1.plot(data1[:,1], yfit(data1[:,1]), label=fitlabel)

logx = np.log(data2[:,2])
logy = np.log(data2[:,1])
coeffs = np.polyfit(logx,logy,deg=1)
poly = np.poly1d(coeffs)
print "Scaling for MPI version: %f" % coeffs[0]

logx = np.log(data3[:,3])
logy = np.log(data3[:,2])
coeffs = np.polyfit(logx,logy,deg=1)
poly = np.poly1d(coeffs)
print "Scaling for Hybrid version: %f" % coeffs[0]

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('System size')
ax1.set_ylabel('Time per step')
ax1.set_title('System size scaling comparison')
ax1.legend(loc='upper left', prop={'size':10})

labels = [''] * len(data1[:,2])
labels[0] = '%d' % int(data1[0,2])
labels[1] = '%d' % int(data1[1,2])
labels[2] = '%d' % int(data1[2,2])
labels[4] = '%d' % int(data1[4,2])
labels[-1] = '%d' % int(data1[-1,2])

ax1.xaxis.set_major_locator(ticker.FixedLocator(data1[:,2]))
ax1.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
ax1.xaxis.set_minor_locator(ticker.NullLocator())

plt.xlim(xmin=min(data1[:,2]), xmax=max(data1[:,2]))

plt.savefig('comparison_system_size_scaling.pdf')
plt.show()
