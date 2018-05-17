import numpy as np
import matplotlib.pyplot as plt


T=np.array(['1'])
init=np.array(['0'])

# T=np.array(['1','1.66','2.33','3'])
# init=np.array(['0','1'])

for index_T, value_T in enumerate(T):
    for index_init, value_init in enumerate(init):
        mat = np.genfromtxt("build/momentaufnahme_" + value_T + "_" + value_init +".mat", unpack=True, skip_header=2)
        plt.pcolormesh(mat, cmap='cool', alpha=1)
        plt.savefig("build/momentaufnahme_" + value_T + "_" + value_init +".pdf")
        plt.close()
        E, m, mm = np.genfromtxt("build/history_" + value_T + "_" + value_init +".mat", unpack=True, skip_header=2)
        plt.plot(E,'o', alpha=0.7, markersize=0.5)
        plt.plot(mm,'o', alpha=0.7, markersize=0.5)
        plt.plot(m,'o', alpha=0.7, markersize=0.5)
        plt.savefig("build/history_" + value_T + "_" + value_init +".pdf")
        plt.close()

T, E_T, m_T, absm_T = np.genfromtxt("build/temperature.mat", unpack=True, skip_header=2)
plt.plot(T, m_T, 'o', label=r'$<m>(T)$', markersize=3, alpha=0.5)
plt.plot(T, absm_T, 'o', label=r'$<|m|>(T)$', markersize=3, alpha=0.5)
#plt.plot(T, E_T, 'o', label=r'E(T)', markersize=3)
plt.axvline(2.27, label='Kritischer Punkt')
plt.xlabel(r'$k_B T$')
plt.legend(loc='best')
plt.savefig('build/temperature.pdf')
