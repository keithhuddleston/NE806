import numpy as np
import matplotlib.pyplot as plt

x, y, = np.loadtxt('ParaStudyNU235.txt')
x2, y2 = np.loadtxt('ParaStudyNU235g2.txt')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 3))
plt.grid(ls='--', axis='y')
ax.plot(x, y, c='r')
plt.ylim(min(y), 1.5)
plt.xlim(min(x), max(x))
plt.xlabel(r'$N({U_{235}})$')
plt.ylabel('Multiplication Factor, k')
ax.plot(x2, y2, c='g')
plt.legend(['Casmo 16', 'Casmo 2'])