import numpy as np
import matplotlib.pyplot as plt
import sys



filename1 = 'strong_scaling/omp_strong_scaling_512.dat'
filename2 = 'strong_scaling/omp_strong_scaling_1024.dat'
filename3 = 'strong_scaling/omp_strong_scaling_2048.dat'
filename4 = 'strong_scaling/omp_strong_scaling_4096.dat'


data1 = np.loadtxt(filename1)
data2 = np.loadtxt(filename2)
data3 = np.loadtxt(filename3)
data4 = np.loadtxt(filename4)

# load data
serial1 = data1[0,1]
serial2 = data2[0,1]
serial3 = data3[0,1]
serial4 = data4[0,1]


if len(data1[0]) > 3: # with error
    plt.errorbar(data1[:,0], serial1/data1[:,1], yerr=data1[:,3], label='N = 512')
    plt.errorbar(data2[:,0], serial2/data2[:,1], yerr=data2[:,3], label='N = 1024')
    plt.errorbar(data3[:,0], serial3/data3[:,1], yerr=data3[:,3], label='N = 2048')
    plt.errorbar(data4[:,0], serial4/data4[:,1], yerr=data4[:,3], label='N = 4096')
else: # without error
    plt.plot(data1[:,0], t1/data1[:,1])#, label='%d mpi-tasks, N = %d' % (data[i*count,0], data[i*count,3]))

# annotate with size
#for i in range(0,len(data)):
#    plt.annotate('%d' % data[i,2], xy=(data[i,0],t1/data[i,1]))


# plot linear speedup line
plt.plot([0,data1[-1,0]+1], [0,data1[-1,0]+1], label='Linear speedup', color='#BDBDBD', linestyle='--')


#size = filename.split('_')[-1].split('.',1)[0]

plt.ylabel(r'Speedup $t_{serial} / t_{parallel}$')
plt.title('Strong scaling of OpenMP version')
plt.xlabel('Number of threads')
plt.legend()
plt.xlim(xmin=0, xmax=data1[-1,0]+1) #, xmax=18)
#plt.ylim(ymax=data[-1,1]+1) #, xmax=18)

plt.savefig('strong_scaling/omp_strong_scaling.pdf')
plt.show()
