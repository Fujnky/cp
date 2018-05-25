import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys

text = sys.stdin.readlines()
q = np.genfromtxt([text[0]], unpack=True, dtype=int)
beta, m, E = np.genfromtxt(text, unpack=True, skip_header=3)

gs = gridspec.GridSpec(2, 1)
ax0 = plt.subplot(gs[0, 0])
ax0.plot(beta, m, 'x', markersize=3)
ax0.set_ylabel(r'$m(\beta)$')
ax1 = plt.subplot(gs[1, 0])
ax1.plot(beta, E, 'x', markersize=3)
ax1.set_ylabel(r'$E(\beta)$')

ax0.axvline(np.log(1+np.sqrt(q)), label='Kritischer Punkt')
ax1.axvline(np.log(1+np.sqrt(q)), label='Kritischer Punkt')


ax1.set_xlabel(r'$\beta$')
ax0.legend(loc='best')
#plt.tight_layout(pad=0.2)
ax0.set_title('q = {}'.format(q))
plt.savefig('build/beta_q{}.pdf'.format(q))
