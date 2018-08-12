import matplotlib.pyplot as plt
import numpy as np
import sys

def plot(file):
    matrix = np.genfromtxt(file)
    theta1, theta2, theta1punkt, theta2punkt = np.split(matrix, 4,  axis=0)
    theta1 = (theta1[0]  + np.pi) % (2 * np.pi) - np.pi
    theta1punkt = theta1punkt[0]
    plt.plot(theta1, theta1punkt, '.', alpha=0.5, markersize=1)
    plt.xlabel(r'$\theta_1$')
    plt.ylabel(r'$\dot\theta_1$')
    plt.tight_layout(pad=0)
    plt.savefig(file + '.pdf')
    plt.close()

plot('build/poincaré_quasi.txt')
plot('build/poincaré_chaos.txt')
plot('build/poincaré_+.txt')


diff_q = np.linalg.norm(np.genfromtxt('build/poincaré_quasi.txt') - np.genfromtxt('build/poincaré_quasi_p.txt'), axis=0)
diff_c = np.linalg.norm(np.genfromtxt('build/poincaré_chaos.txt') - np.genfromtxt('build/poincaré_chaos_p.txt'), axis=0)
t = np.array(range(len(diff_q)))/50
plt.plot(t, diff_q, '.', markersize=1, alpha=0.5, label='Quasiperiodisch')
plt.plot(t, diff_c, '.', markersize=1, alpha=0.5, label='Chaotisch')
plt.ylabel('Euklid. Abst. des jew. Phasenraumpunkts zum ungestörten Syst.')
plt.xlabel('$t$/s')
plt.legend(loc='best')
plt.yscale('log')
plt.tight_layout(pad=0)
plt.savefig('build/euclid.pdf')
