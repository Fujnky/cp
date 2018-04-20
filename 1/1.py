import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

N, rN = np.genfromtxt('build/output.txt', unpack=True)
rN = np.sqrt(rN)

lr = LinearRegression()
lr.fit(np.log(N.reshape(-1, 1)), np.log(rN))

plt.plot(N, rN, 'C0x')
N_ = np.linspace(10, 60, 1000)
print(lr.coef_)
plt.plot(N_, N_**(lr.coef_[0])*np.exp(lr.intercept_), 'C1-')
#plt.ylim(0, 100)
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel("$N$")
plt.ylabel("$R_N$")
plt.savefig('build/1.pdf')
