import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys


text = sys.stdin.readlines()
sizex, sizey, sweeps = np.genfromtxt([text[1]], unpack=True, dtype=int)
data = np.genfromtxt(text, unpack=True, skip_header=2).reshape(sizex, sweeps, sizey)

time = sweeps-1

fig  = plt.figure()
ax   = plt.subplot(111)
mesh = ax.pcolormesh(data[:,0,:], cmap='viridis')

def animate(t):
    print('\r{}'.format(t), end='')
    ax.cla()
    mesh = ax.pcolormesh(data[:,t+1,:], cmap='viridis')
    ax.set_title('Wolff – Step: {}'.format(t*10), loc='right')
    return mesh

#writer = animation.ImageMagickWriter(fps=25)
ani = animation.FuncAnimation(fig, animate, frames=time, interval=0)
#ani.save('build/anim.mp4')
plt.show()
