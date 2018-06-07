import matplotlib.pyplot as plt
import numpy as np

matrix = np.genfromtxt('build/hamiltonian', unpack=True)
plt.pcolormesh(matrix)
plt.gca().invert_yaxis()
plt.colorbar()
plt.show()
