import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def Energieerhaltung(r, v):
    for index, value in enumerate(r[0, :]):
        E_kin[index] = 0.5 * v[index]
        E_pot[index] = 0.5 * r[index]
    return E_kin+E_pot


matrix = np.genfromtxt('build/rk4_kepler.txt')

r, v = np.split(matrix, 2,  axis=0)

E_kin = 0.5 * np.linalg.norm(v, axis=0)**2
E_pot = 0.5 * np.linalg.norm(r, axis=0)**2
E_pot_kepler = 1 / np.linalg.norm(r, axis=0)


#plt.plot(E_kin)
#plt.plot(E_pot_kepler)
plt.plot(E_kin + E_pot_kepler)
# plt.show()
plt.close()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(*r)
# plt.show()
plt.close()

r_shifted = r[:,1:]
r = r[:,:-1]
a = np.linalg.norm(r, axis=0)
b = np.linalg.norm(r_shifted, axis=0)
c = np.linalg.norm(r - r_shifted, axis=0)
area = np.sqrt(4 * a**2 * b**2 - (a**2 + b**2 - c**2)**2)/4
plt.plot(1 - area/area.mean())
plt.show()
