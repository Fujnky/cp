import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def funktion(N,a,x):
    return a*N**x


matrix = np.genfromtxt('build/submatrix', unpack=True)
plt.pcolormesh(matrix, cmap='seismic')
plt.gca().invert_yaxis()
plt.colorbar()
plt.savefig('build/matrix.pdf')
plt.close()


ev = np.genfromtxt('build/ev', unpack=True)
ev.sort()
plt.plot(ev, 'x', markersize=2, alpha=0.25)
plt.savefig('build/eigenvalues.pdf')
plt.close()

N, t = np.genfromtxt('build/runtimes', unpack=True)
N_ = np.linspace(2, 14, 100)
params, pcov = curve_fit(funktion, N, t, maxfev=900000000)
plt.plot(N_, funktion(N_, *params, ))
plt.plot(N, t, 'x')
plt.title("x="+str(params[1]))
# plt.yscale('log')
plt.savefig('build/runtime.pdf')
