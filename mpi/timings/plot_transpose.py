import numpy as np
import matplotlib.pyplot as plt
import sys

# get filename
#if len(sys.argv) < 2:
#       filename = raw_input('Data-file name: ')
#else:
#       filename = sys.argv[1]


size = 1024

filename1 = 'daint_mpi_datatype_%d.dat' % size
filename2 = 'daint_mpi_localtranspose_%d.dat' % size
#filename2 = 'daint_mpi_datatype.dat'


# load data
data1 = np.loadtxt(filename1) # number of cpus, time, size
data2 = np.loadtxt(filename2) # number of cpus, time, size

plt.plot(data1[:,0], data1[:,1], '-o', label='Datatype')
plt.plot(data2[:,0], data2[:,1], '-o', label='Local transpose')

# annotate with size
for i in range(0,len(data1)):
        if data1[i,0] == data2[i,0]:
                plt.annotate('%d' % data1[i,2], xy=(data1[i,0], (data1[i,1]+data2[i,1])/2))
        else:
                plt.annotate('%d' % data1[i,2], xy=(data1[i,0],data1[i,1]))
#for i in range(0,len(data2)):
#       plt.annotate('%d' % data2[i,2], xy=(data2[i,0],data2[i,1]))


plt.xlabel('number of CPUs')
plt.ylabel('time')
plt.title('MPI on Piz Daint')
plt.legend()
plt.xlim(xmin=0, xmax=max(data1[-1,0],data2[-1,0])+1)


#size = filename1.split('_')[-1].split('.',1)[0]

plt.savefig('daint_mpi_transpose_%d.pdf' % size)
plt.show()
