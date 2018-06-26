import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import sys

text = sys.stdin.readlines()
matrix = np.genfromtxt(text)
theta2, theta1, theta2punkt, theta1punkt = np.split(matrix, 4,  axis=0)

theta1 = theta1[0]
theta2 = theta2[0]
theta1punkt = theta1punkt[0]
theta2punkt = theta2punkt[0]

plt.plot(theta1)
plt.plot(theta2)

plt.tight_layout(pad=0)
plt.savefig('build/1.pdf')
plt.close()

m1 = 1
m2 = 1

L1 = 1
L2 = 1

g = 9.81

T = 0.5 * m1 * theta1punkt**2 * L1**2 + 0.5* (m2 * theta2punkt**2 * L2**2 + m2 * theta1punkt**2 * L1**2 + 2* m2 * theta1punkt * L1 * theta2punkt * L2 * np.cos(theta1-theta2))

V = -m1 * g * L1 * np.cos(theta1) - m2 * g * L1 * np.cos(theta1) - m2 * g * L2 * np.cos(theta2)

plt.plot(T+V, label='T+V')
#plt.plot(theta2)
plt.legend(loc='best')
plt.tight_layout(pad=0)
plt.savefig('build/energie.pdf')
plt.close()

fig  = plt.figure()
ax   = plt.subplot(111)

def animate(t):
    print('\r{}'.format(t), end='')
    ax.cla()
    plot = ax.plot((0, L1*np.sin(theta1[t]), L1*np.sin(theta1[t]) + L2*np.sin(theta2[t])), (0, -L1*np.cos(theta1[t]), -L1*np.cos(theta1[t])-L2*np.cos(theta2[t])))
    plt.axis('equal')
    plt.axis([-4, 4, -4, 4])
    ax.set_title('Doppelpendel – Step: {}'.format(t), loc='right')
    return plot

#writer = animation.ImageMagickWriter(fps=25)
ani = animation.FuncAnimation(fig, animate, frames=theta1.shape[0])
#plt.show()
#ani.save('build/anim.mp4', fps=50)
#print()
