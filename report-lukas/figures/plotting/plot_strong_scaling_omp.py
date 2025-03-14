import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import sys



# load data
filename1 = 'omp_strong_scaling_512.dat'
filename2 = 'omp_strong_scaling_1024.dat'
filename3 = 'omp_strong_scaling_2048.dat'
filename4 = 'omp_strong_scaling_4096.dat'


data1 = np.loadtxt(filename1)
data2 = np.loadtxt(filename2)
data3 = np.loadtxt(filename3)
data4 = np.loadtxt(filename4)

# not actually serial data, but omp version with one thread
serial1 = data1[0,1]
serial2 = data2[0,1]
serial3 = data3[0,1]
serial4 = data4[0,1]

fig1, ax1 = plt.subplots()

# plot
if len(data1[0]) > 3: # with error
    ax1.errorbar(data1[:,0], serial1/data1[:,1], yerr=data1[:,3], label='N = 512', color='r')
    ax1.errorbar(data2[:,0], serial2/data2[:,1], yerr=data2[:,3], label='N = 1024', color='b')
    ax1.errorbar(data3[:,0], serial3/data3[:,1], yerr=data3[:,3], label='N = 2048', color='g')
    ax1.errorbar(data4[:,0], serial4/data4[:,1], yerr=data4[:,3], label='N = 4096', color='y')
else: # without error
    ax1.plot(data1[:,0], t1/data1[:,1])#, label='%d mpi-tasks, N = %d' % (data[i*count,0], data[i*count,3]))

# annotate with size
#for i in range(0,len(data)):
#    ax1.annotate('%d' % data[i,2], xy=(data[i,0],t1/data[i,1]))


# plot linear speedup line
ax1.plot([0,data1[-1,0]+1], [0,data1[-1,0]+1], label='Linear speedup', color='k', linestyle='--')


#size = filename.split('_')[-1].split('.',1)[0]

ax1.set_ylabel(r'Speedup $t_1 / t_n$')
ax1.set_title('OpenMP version strong scaling')
ax1.set_xlabel('Number of threads')
ax1.legend()
ax1.set_xlim(xmin=0, xmax=data1[-1,0]+1) #, xmax=18)
#plt.ylim(ymax=data[-1,1]+1) #, xmax=18)

labels = ['%d' % d for d in data1[:,0]]
#labels = int(data1[:,0])
#labels[0] = '%d' % int(data1[0,0])
#labels[-4] = '%d' % int(data1[-4,0])
#labels[-3] = '%d' % int(data1[-3,0])
#labels[-2] = '%d' % int(data1[-2,0])
#labels[-1] = '%d' % int(data1[-1,0])

ax1.xaxis.set_major_locator(ticker.FixedLocator(data1[:,0]))
ax1.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
ax1.xaxis.set_minor_locator(ticker.NullLocator())

plt.savefig('omp_strong_scaling.pdf')
plt.show()
