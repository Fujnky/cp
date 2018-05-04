import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def plot(idx, min_, max_, func):
    numbers = np.genfromtxt('build/output_'+idx+'.txt')
    plt.hist(numbers, bins=np.linspace(min_, max_, 30), density=True)
    x = np.linspace(min_, max_, 1000)
    plt.plot(x, func(x))
    plt.xlim(min_, max_)
    plt.savefig('build/a2_hist_' + idx + '.pdf')
    plt.close()

def gauss(x,sigma,mu):
    return (1/np.sqrt(2*np.pi*sigma**2))*np.exp(-(x-mu)**2/(2*sigma**2))

plot("1", 0, np.pi/2, np.cos)
plot("2", 0, np.pi/2, np.cos)
plot("3", -7, 13, lambda x: gauss(x, 2, 3))
