import numpy as np
import matplotlib.pyplot as plt


def mag(H):
    return np.tanh(H)

#laden
H, m = np.genfromtxt('build/ising0d.txt', unpack=True)

plt.plot(H, m, 'C0.', alpha=0.3, label='Simulation')
plt.plot(H,mag(H),'C1-', label=r'$\tanh(\beta H)$')

plt.xlabel("H")
plt.ylabel("m")
plt.legend(loc='best')
plt.savefig('build/1.pdf')
