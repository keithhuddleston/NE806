import numpy as np
import matplotlib.pyplot as plt

x, y, = np.loadtxt('ParaStudyNU235.txt')


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 3))
plt.grid(ls='--', axis='y')
ax.plot(x, y, c='k')
plt.ylim(min(y), 1.05)
plt.xlim(min(x), max(x))
plt.xlabel(r'$N({U_{235}})$')
plt.ylabel('Multiplication Factor, k')