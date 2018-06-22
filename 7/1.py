import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

matrix = np.genfromtxt('build/rk4.txt')
matrix2 = np.genfromtxt('build/rk4_har.txt')
matrix3 = np.genfromtxt('build/rk4_np.txt')


r, v = np.split(matrix, 2,  axis=0)
r2, v2 = np.split(matrix2, 2,  axis=0)
r3, v3 = np.split(matrix3, 2,  axis=0)


E_kin = 0.5 * np.linalg.norm(v, axis=0)**2
E_pot = 0.5 * np.linalg.norm(r, axis=0)**2

plt.plot(E_kin, label=r'$E_\mathrm{kin}$')
plt.plot(E_pot, label=r'$E_\mathrm{pot}$')
plt.plot(E_kin + E_pot, label=r'$E_\mathrm{kin} + E_\mathrm{pot}$')
plt.legend(loc='best')
plt.tight_layout(pad=0)
plt.savefig('build/energie.pdf')
plt.close()


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(*r , label='a)')
ax.plot(*r2, label='b) 1')
ax.plot(*r3, label='b) 2')
plt.legend(loc='best')
plt.tight_layout(pad=0)
plt.savefig('build/plot.pdf')
plt.close()
