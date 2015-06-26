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
data1 = np.loadtxt('hybrid/hybrid_size_data.dat') # ntasks, nthreads, time, size, err
data2 = np.loadtxt('hybrid/hybrid_size_lt.dat') # ntasks, nthreads, time, size, err

fig1, ax1 = plt.subplots()

if len(data1[0]) > 2: # with error
    ax1.errorbar(data1[:,3], data1[:,2], yerr=data1[:,4], label='Transpose with data type')
    ax1.errorbar(data2[:,3], data2[:,2], yerr=data2[:,4], label='Locally transpose blocks')
else: # without error
    ax1.plot(data1[:,1], data1[:,0], label='vectorized (-ftree-vectorize)')

# annotate with size
#for i in range(0,len(data)):
#    plt.annotate('%d' % data[i,2], xy=(data[i,0],t1/data[i,1]))

#size = filename.split('_')[-1].split('.',1)[0]

logx = np.log(data1[:,3])
logy = np.log(data1[:,2])
coeffs = np.polyfit(logx,logy,deg=1)
poly = np.poly1d(coeffs)

yfit = lambda x: np.exp(poly(np.log(x)))

#print "Scaling for serial version: %f" % coeffs[0]
#fitlabel = r'$\alpha =$%f' % coeffs[0]
#ax1.plot(data1[:,1], yfit(data1[:,1]), label=fitlabel)


ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('System size')
ax1.set_ylabel('Time per step')
ax1.set_title('System size scaling of hybrid version\n(32 MPI tasks, 8 OpenMP threads)')
ax1.legend(loc='upper left')

labels = [''] * len(data1[:,3])
labels[0] = '%d' % int(data1[0,3])
labels[1] = '%d' % int(data1[1,3])
labels[2] = '%d' % int(data1[2,3])
labels[4] = '%d' % int(data1[4,3])
labels[-1] = '%d' % int(data1[-1,3])

ax1.xaxis.set_major_locator(ticker.FixedLocator(data1[:,3]))
ax1.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
ax1.xaxis.set_minor_locator(ticker.NullLocator())

plt.xlim(xmin=min(data1[:,3]), xmax=max(data1[:,3]))

plt.savefig('hybrid/hybrid_size_scaling.pdf')
plt.show()
