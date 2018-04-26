import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def plot(filename, period):
    numbers = np.genfromtxt(filename+'.txt')
    plt.hist(numbers, bins=10)
    plt.savefig(filename+'_hist.pdf')
    plt.close()

    plt.plot(numbers[:min(period, numbers.shape[0])-1], numbers[1:min(period, numbers.shape[0])], '.', markersize=2)
    plt.savefig(filename+'_scatter.pdf')
    plt.close()

#plot('build/1', 6075)
#plot('build/2', 256)
#plot('build/3', 2147483648)
#plot('build/4', 2147483647)
#plot('build/5', 100000)
#plot('build/6', 100000)

data = np.genfromtxt('build/period.txt', unpack=True)
x = np.arange(0.5, 16.5)
plt.pcolormesh(x, x, data, norm=colors.LogNorm(vmin=data.min(), vmax=data.max()))
plt.xlabel('b')
plt.ylabel('c')
plt.colorbar(label='period')
plt.tight_layout(pad=0)
plt.savefig('build/period.pdf')
plt.close()
