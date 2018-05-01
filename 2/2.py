import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def plot(idx):
    numbers = np.genfromtxt('build/output_'+idx+'.txt')
    plt.hist(numbers, bins=10)
    plt.savefig('build/a2_hist_' + idx + '.pdf')
    plt.close()

plot("1")
plot("2")
