import numpy as np
import matplotlib.pyplot as plt

mat = np.genfromtxt('build/momentaufnahme.mat', unpack=True, skip_header=2)

plt.pcolormesh(mat, cmap='coolwarm')
plt.savefig('build/momentaufnahme1.pdf')
