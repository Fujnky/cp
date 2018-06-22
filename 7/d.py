import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys

text = sys.stdin.readlines()
matrix = np.genfromtxt(text)
r1,r2,r3,v1,v2,v3 = np.split(matrix, 6,  axis=0)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(*r1)
ax.plot(*r2)
ax.plot(*r3)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.tight_layout(pad=0)
plt.savefig('build/threebody.pdf')
plt.close()
