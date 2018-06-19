import numpy as np
import matplotlib.pyplot as plt

lanczos = np.genfromtxt('build/lanczos.txt', unpack=True)
householder = np.genfromtxt('build/householder.txt', unpack=True)
eigen = np.genfromtxt('build/eigen.txt', unpack=True)
N = np.linspace(3, 20, 18)

plt.plot(N, lanczos, 'x', label='Lanczos')
plt.plot(N, householder, 'x', label='Householder')
plt.plot(N, eigen, 'x', label='Eigen')
plt.xlabel('N')
plt.ylabel('t / µs')
plt.legend(loc='best')
plt.tight_layout(pad=0)
plt.savefig('build/1.pdf')
