import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys


text = sys.stdin.readlines()
save_every_nth_step, q = np.genfromtxt([text[0]], unpack=True, dtype=int)
sizex, sizey, steps = np.genfromtxt([text[2]], unpack=True, dtype=int)
data = np.genfromtxt(text, unpack=True, skip_header=3).reshape(sizex, steps, sizey)

time = steps-1

fig  = plt.figure()
ax   = plt.subplot(111)
mesh = ax.pcolormesh(data[:,0,:], cmap='viridis')

def animate(t):
    print('\r{}'.format(t), end='')
    ax.cla()
    mesh = ax.pcolormesh(data[:,t+1,:], cmap='viridis')
    ax.set_title('Wolff-Algorithmus für q = {} – Step: {}'.format(q, t*save_every_nth_step), loc='right')
    if t == 0 or t == int(time/2) or t == time-1:
        plt.savefig('build/q{}_t{}.pdf'.format(q, int(2*t/(time-1))))
    return mesh

#writer = animation.ImageMagickWriter(fps=25)
ani = animation.FuncAnimation(fig, animate, frames=time)
ani.save('build/anim_q{}.mp4'.format(q), fps=15)
print()
#plt.show()
