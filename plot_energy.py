import matplotlib.pyplot as plt
import matplotlib.ticker as ticker



cpus = [1, 2, 4, 8, 16, 32, 64, 128, 256]
energy = [66982, 30434, 36037, 25470, 34324, 13861, 17199, 32931, 52283]

fig1, ax1 = plt.subplots()

ax1.plot(cpus, energy)

ax1.set_title("Energy consumption of MPI version")
ax1.set_xlabel("Number of CPU's")
ax1.set_ylabel("Energy consumption in joules")

labels = [''] * len(cpus)
labels[-1] = '%d' % int(cpus[-1])
labels[-2] = '%d' % int(cpus[-2])
labels[-3] = '%d' % int(cpus[-3])
labels[-4] = '%d' % int(cpus[-4])
labels[-5] = '%d' % int(cpus[-5])

ax1.xaxis.set_major_locator(ticker.FixedLocator(cpus))
ax1.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
ax1.xaxis.set_minor_locator(ticker.NullLocator())


plt.savefig("energy.pdf")
plt.show()

