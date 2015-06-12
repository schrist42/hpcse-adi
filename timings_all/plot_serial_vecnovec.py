import numpy as np
import matplotlib.pyplot as plt
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


if len(data_vec[0]) > 2: # with error
    plt.errorbar(data_vec[:,1], data_vec[:,0], yerr=data_vec[:,2], label='vectorized (-ftree-vectorize)')
    plt.errorbar(data_novec[:,1], data_novec[:,0], yerr=data_novec[:,2], label='not vectorized (-fno-tree-vectorize)')
else: # without error
    plt.plot(data_vec[:,1], data_vec[:,0], label='vectorized (-ftree-vectorize)')
    plt.plot(data_novec[:,1], data_novec[:,0], label='not vectorized (-fno-tree-vectorize)')

# annotate with size
#for i in range(0,len(data)):
#    plt.annotate('%d' % data[i,2], xy=(data[i,0],t1/data[i,1]))

#size = filename.split('_')[-1].split('.',1)[0]

plt.xlabel('size')
plt.ylabel('time per step')
plt.title('Comparison vectorizing and not vectorizing')
plt.legend()
#plt.xlim(xmin=0, xmax=data_vec[-1,0]+1) #, xmax=18)

plt.savefig('serial_vecnovec.pdf')
plt.show()
