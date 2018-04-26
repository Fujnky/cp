import numpy as np
import matplotlib.pyplot as plt

def plot(filename, period):
    numbers = np.genfromtxt(filename+'.txt')
    plt.hist(numbers, bins=10)
    plt.savefig(filename+'_hist.pdf')
    plt.close()

    plt.plot(numbers[:min(period, numbers.shape[0])-1], numbers[1:min(period, numbers.shape[0])], '.', markersize=2)
    plt.savefig(filename+'_scatter.pdf')
    plt.close()

plot('build/1', 6075)
plot('build/2', 256)
plot('build/3', 2147483648)
plot('build/4', 2147483647)
