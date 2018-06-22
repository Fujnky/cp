import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

matrix = np.genfromtxt('build/rk4_kepler.txt')
r, v = np.split(matrix, 2,  axis=0)

matrix = np.genfromtxt('build/rk4_kepler2.txt')
r2, v2 = np.split(matrix, 2,  axis=0)

matrix = np.genfromtxt('build/rk4_kepler_alpha_0.9.txt')
r3, v3 = np.split(matrix, 2,  axis=0)

matrix = np.genfromtxt('build/rk4_kepler_alpha_1.1.txt')
r4, v4 = np.split(matrix, 2,  axis=0)

E_kin = 0.5 * np.linalg.norm(v, axis=0)**2
E_pot_kepler = -1 / np.linalg.norm(r, axis=0)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(*r, label=r'$\alpha = 1$')
#shift for visibility
r3[2,:] += 1
r4[2,:] -= 1
ax.plot(*r3, label=r'$\alpha = 0.9$')
ax.plot(*r4, label=r'$\alpha = 1.1$')
plt.legend(loc='best')
plt.tight_layout(pad=0)
plt.savefig('build/kepler_plot.pdf')
plt.close()


plt.plot(E_kin, label=r'$E_\mathrm{kin}$')
plt.plot(E_pot_kepler, label=r'$V$')
plt.plot(E_kin + E_pot_kepler, label=r'$E_\mathrm{kin} + V$')
plt.legend(loc='best')
plt.tight_layout(pad=0)
plt.savefig('build/kepler_energie.pdf')
plt.close()


def calc(r, v):
    a = np.linalg.norm(r.T - r[:,0], axis=1)
    revs = np.argwhere(np.r_[True, a[1:] < a[:-1]] & np.r_[a[:-1] < a[1:], True])
    T = np.diff(revs[:-1,0]).mean()
    r_ = np.linalg.norm(r, axis=0)
    return T**2 / (r_.max() + r_.min())**3

print('Gleich wenn Kepler 3 erfüllt:')
print(calc(r, v), calc(r2, v2))
print('Juhu')

r_shifted = r[:,1:]
r = r[:,:-1]
a = np.linalg.norm(r, axis=0)
b = np.linalg.norm(r_shifted, axis=0)
c = np.linalg.norm(r - r_shifted, axis=0)
area = np.sqrt(4 * a**2 * b**2 - (a**2 + b**2 - c**2)**2)/4
plt.plot(1 - area/area.mean(), label='Abweichung der Dreiecksfläche,\ndie pro Zeiteinheit überstrichen\nwird, vom Mittelwert.')
plt.legend(loc='best')
plt.tight_layout(pad=0)
plt.savefig('build/kepler2.pdf')
plt.close()
