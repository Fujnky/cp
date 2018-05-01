import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

#laden
N, rN = np.genfromtxt('build/output_.txt', unpack=True)
N2, rN2 = np.genfromtxt('build/output.txt', unpack=True)
#gespeichert ist R_N^2
rN = np.sqrt(rN)
rN2 = np.sqrt(rN2)

#Lineare Regression auf den logarithmierten Daten um den Exponenten zu ermitteln
lr = LinearRegression()
lr.fit(np.log(N.reshape(-1, 1)), np.log(rN))

plt.plot(N, rN, 'C0x')
plt.plot(N2, rN2, 'C2x')

N_ = np.linspace(5, 1450, 1000)
print(lr.coef_)
plt.plot(N_, N_**(lr.coef_[0])*np.exp(lr.intercept_), 'C1-')
#plt.ylim(0, 100)
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel("$N^{{{0:3.2f}}}$".format(lr.coef_[0]))
plt.ylabel("$R_N$")
plt.savefig('build/1.pdf')
