import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.rcParams["figure.figsize"] = (7,8)


T=np.array(['1.00', '1.66', '2.33', '3.00'])
init=np.array(['0', '1'])

# T=np.array(['1','1.66','2.33','3'])
# init=np.array(['0','1'])

def multiplot(name, temps, states, E, m, absm):
    gs = gridspec.GridSpec(len(temps), 3)
    for i, temp in enumerate(temps):
        ax0 = plt.subplot(gs[i, 0])
        ax0.pcolormesh(states[i], cmap='cool', alpha=1)
        ax0.set_ylabel(r'$k_BT = {}$'.format(temp))
        plt.setp(ax0.get_xticklabels(), visible=False)
        plt.setp(ax0.get_yticklabels(), visible=False)

        ax1 = plt.subplot(gs[i, 1])
        ax1.plot(E[i],'x', alpha=0.7, markersize=0.7)

        ax2 = plt.subplot(gs[i, 2])
        ax2.plot(m[i],'o', alpha=0.7, markersize=0.7, label='$<m>$')
        ax2.plot(absm[i],'o', alpha=0.7, markersize=0.7, label='$<|m|>$')

        if i == 0:
            ax2.legend(loc='best')
            ax0.set_title('Momentaufnahme')
            ax1.set_title('Energie')
            ax2.set_title('Magnetisierung')

        if i != len(temps)-1:
            plt.setp(ax1.get_xticklabels(), visible=False)
            plt.setp(ax2.get_xticklabels(), visible=False)
        else:
            ax1.set_xlabel('t')
            ax2.set_xlabel('t')


    plt.tight_layout(pad=0)
    plt.savefig('build/multiplot{}.pdf'.format(name))
    plt.close()




for index_init, value_init in enumerate(init):
    mat = np.zeros((len(T), 100, 100))
    E = np.zeros((len(T), 1000))
    m = absm = np.zeros((len(T), 1000))
    for index_T, value_T in enumerate(T):
        mat[index_T] = np.genfromtxt("build/momentaufnahme_" + value_T + "_" + value_init +".mat", unpack=True, skip_header=2)
        E[index_T], m[index_T], absm[index_T] = np.genfromtxt("build/history_" + value_T + "_" + value_init +".mat", unpack=True, skip_header=2)
    multiplot(value_init, T, mat, E, m, absm)


##Hier alles von T abhängige

T, E_T, m_T, absm_T, c_T = np.genfromtxt("build/temperature.mat", unpack=True, skip_header=2)
gs = gridspec.GridSpec(4, 1)
ax0 = plt.subplot(gs[0, 0])
ax0.plot(T, m_T, 'x', markersize=3)
ax0.set_ylabel('$<m>(T)$')
ax1 = plt.subplot(gs[1, 0])
ax1.plot(T, absm_T, 'x', markersize=3)
ax1.set_ylabel('$<|m|>(T)$')
ax2 = plt.subplot(gs[2, 0], sharex=ax0)
ax2.plot(T, E_T, 'x', markersize=3)
ax2.set_ylabel('$E(T)$')
ax3 = plt.subplot(gs[3, 0], sharex=ax0)
ax3.plot(T, c_T, 'x', markersize=3)
ax3.set_ylabel('$c(T)$')
plt.setp(ax0.get_xticklabels(), visible=False)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
ax0.axvline(2.27, label='Kritischer Punkt')
ax1.axvline(2.27, label='Kritischer Punkt')
ax2.axvline(2.27, label='Kritischer Punkt')
ax3.axvline(2.27, label='Kritischer Punkt')
ax3.set_xlabel(r'$k_B T$')
ax0.legend(loc='best')
plt.tight_layout(pad=0)
plt.savefig('build/temperature.pdf')
