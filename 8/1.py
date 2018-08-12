import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
from matplotlib.collections import LineCollection
import matplotlib.colors
import matplotlib.patches as mpatches
import numpy as np
import sys

steps_per_second = 50
t_max = 40
steps = steps_per_second * t_max

matrix = np.genfromtxt('build/eins.txt')
theta1, theta2, theta1punkt, theta2punkt = np.split(matrix, 4,  axis=0)

#frag mich nicht
theta1 = theta1[0]
theta2 = theta2[0]
theta1punkt = theta1punkt[0]
theta2punkt = theta2punkt[0]

m1 = 1
m2 = 1

L1 = 0.4
L2 = 1

g = 9.81

#energien ausrechnen
T1 = 0.5 * m1 * theta1punkt**2 * L1**2
T2 = 0.5* (m2 * theta2punkt**2 * L2**2 + m2 * theta1punkt**2 * L1**2 + 2* m2 * theta1punkt * L1 * theta2punkt * L2 * np.cos(theta1-theta2))
V1 = -m1 * g * L1 * np.cos(theta1)
V2 = - m2 * g * L1 * np.cos(theta1) - m2 * g * L2 * np.cos(theta2)

#pot umeichen für barplots
gauge = np.min((V1, V2))
V1g = V1-gauge
V2g = V2-gauge
max_e = np.max((T1, T2, V1, V2))


fig  = plt.figure()
gs = gridspec.GridSpec(2, 2, height_ratios=[4, 1], width_ratios=[4, 1])
ax1 = plt.subplot(gs[0, 0])
ax2 = plt.subplot(gs[1, :])
ax3 = plt.subplot(gs[0, 1])

#spaß mit linien
x = L1*np.sin(theta1) + L2*np.sin(theta2)
y = -L1*np.cos(theta1)-L2*np.cos(theta2)
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
cm = matplotlib.colors.LinearSegmentedColormap.from_list('own', [(255/255, 127/255, 14/255, 0), (255/255, 127/255, 14/255, 1)])
lc = LineCollection(segments, cmap=cm)
lc.set_linewidth(2)
ax1.add_collection(lc)

#pendel
plot, = ax1.plot(0, '.-C0', markersize=10)
ax1.axis('equal')
ax1.axis([-3, 3, -3, 3])

#energieverlauf
plot_T, = ax2.plot(0, 'C0', label=r'$T$')
plot_V, = ax2.plot(0, 'C1', label=r'$V$')
plot_TV, = ax2.plot(0, 'C2', label='T+V')
ax2.set_xlim(0, t_max)
all_values = np.concatenate((T1+T2, V1+V2, T1+V1+V2+T2))
ax2.set_ylim(all_values.min()-5, all_values.max()+5)
ax2.legend(loc='center right')
ax2.set_ylabel("$E$/J")
ax2.set_xlabel("$t$/s")

#barplots
bars = ax3.bar([0.8,1.2,1.8,2.2], (0,0,0,0), 0.4)
bars[1].set_color('C1')
bars[3].set_color('C1')
ax3.set_xticklabels((0, '$m_1$', '$m_2$'))
ax3.set_ylim(0, max_e*1.2)
ax3.legend((bars[0], bars[1]), (r'$E_\mathrm{kin}$', r'$E_\mathrm{pot}$'))
ax3.yaxis.set_label_position("right")
ax3.yaxis.tick_right()
ax3.set_ylabel("$E$/J")

#anmimieren
def animate(t):
    print('\r{}'.format(t), end='')
    color_array = np.zeros(t_max * steps_per_second)
    if t > 2: #frag mich nicht
        #neue farben berechnen für schweif
        num = int(1.5*steps_per_second) if t > int(1.5*steps_per_second) else t
        color_array[t-num:t-1] = np.linspace(0, 1, num-1)
        lc.set_array(color_array)

    #pendel setzen
    x, y = ((0, L1*np.sin(theta1[t-1]), L1*np.sin(theta1[t-1]) + L2*np.sin(theta2[t-1])),
    (0, -L1*np.cos(theta1[t-1]), -L1*np.cos(theta1[t-1])-L2*np.cos(theta2[t-1])))
    plot.set_xdata(x)
    plot.set_ydata(y)

    #energieverlauf setzen
    T1_ = T1[:t]
    T2_ = T2[:t]
    V1_ = V1[:t]
    V2_ = V2[:t]
    plot_T.set_xdata(np.array(range(len(V1_)))/steps_per_second)
    plot_T.set_ydata(T1_+T2_)
    plot_V.set_xdata(np.array(range(len(V1_)))/steps_per_second)
    plot_V.set_ydata(V1_+V2_)
    plot_TV.set_xdata(np.array(range(len(V1_)))/steps_per_second)
    plot_TV.set_ydata(T1_+V1_+T2_+V2_)

    #bars setzen
    bars[0].set_height(T1[t])
    bars[1].set_height(V1g[t])
    bars[2].set_height(T2[t])
    bars[3].set_height(V2g[t])

    return [plot, plot_T, plot_V, plot_TV, lc]+[bar for bar in bars]

def init():
    return animate(0)

ani = animation.FuncAnimation(fig, animate, frames=range(3, theta1.shape[0]), interval=1.0/steps_per_second*1000, init_func=init, blit=True)
plt.tight_layout(pad=.5)
plt.show()
#ani.save('build/anim.mp4', fps=steps_per_second, bitrate=-1, dpi=225)
