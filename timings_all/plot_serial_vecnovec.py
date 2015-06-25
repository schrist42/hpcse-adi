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
data_vec = np.loadtxt('serial/serial_vec.dat') # time, size, error
data_novec = np.loadtxt('serial/serial_novec.dat') # time, size, error

fig1, ax1 = plt.subplots()

if len(data_vec[0]) > 2: # with error
    ax1.errorbar(data_vec[:,1], data_vec[:,0], yerr=data_vec[:,2], label=r'vectorized (-ftree-vectorize)')
    ax1.errorbar(data_novec[:,1], data_novec[:,0], yerr=data_novec[:,2], label=r'not vectorized (-fno-tree-vectorize)')
else: # without error
    ax1.plot(data_vec[:,1], data_vec[:,0], label='vectorized (-ftree-vectorize)')
    ax1.plot(data_novec[:,1], data_novec[:,0], label='not vectorized (-fno-tree-vectorize)')

# annotate with size
#for i in range(0,len(data)):
#    plt.annotate('%d' % data[i,2], xy=(data[i,0],t1/data[i,1]))

#size = filename.split('_')[-1].split('.',1)[0]

logx = np.log(data_vec[:,1])
logy = np.log(data_vec[:,0])
coeffs = np.polyfit(logx,logy,deg=1)
poly = np.poly1d(coeffs)

yfit = lambda x: np.exp(poly(np.log(x)))

print "Scaling for serial version: %f" % coeffs[0]
fitlabel = r'$\alpha =$%f' % coeffs[0]
ax1.plot(data_vec[:,1], yfit(data_vec[:,1]), label=fitlabel)


ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('Size')
ax1.set_ylabel('Time per step')
ax1.set_title('System size scaling of serial version')
ax1.legend(loc='upper left')

labels = [''] * len(data_vec[:,1])
labels[0] = '%d' % int(data_vec[0,1])
labels[1] = '%d' % int(data_vec[1,1])
labels[2] = '%d' % int(data_vec[2,1])
labels[4] = '%d' % int(data_vec[4,1])
labels[-1] = '%d' % int(data_vec[-1,1])

ax1.xaxis.set_major_locator(ticker.FixedLocator(data_vec[:,1]))
ax1.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
ax1.xaxis.set_minor_locator(ticker.NullLocator())

plt.xlim(xmin=min(data_vec[:,1]), xmax=max(data_vec[:,1]))

plt.savefig('serial/serial_vecnovec.pdf')
plt.show()
