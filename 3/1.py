import numpy as np
import matplotlib.pyplot as plt


def mag(H):
    return np.tanh(H)

#laden
H, m = np.genfromtxt('build/ising0d.txt', unpack=True)

plt.plot(H, m, 'C0.', alpha=0.2)
plt.plot(H,mag(H),'C1-')

plt.xlabel("H")
plt.ylabel("m")
plt.savefig('build/1.pdf')
