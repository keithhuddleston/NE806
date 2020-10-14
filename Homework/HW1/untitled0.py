from scipy.special import expn
import numpy as np
import matplotlib.pyplot as plt

source = 2
sigmaa = 5
width = 1


def psiplus1(x, mu, source, sigmaa):
    return source/(2*sigmaa) * (1 - np.exp(-sigmaa*x/mu))


def psiplus2(x, mu, source, sigmaa, width):
    return source/(2*sigmaa) * (-1 + np.exp(sigmaa*width/(2*mu))) * np.exp(-sigmaa*x/mu)


def psimin(x, mu, source, sigmaa, width):
    return source/(2*sigmaa) * (1 - np.exp(-sigmaa*(width/2-x)/mu))


x1 = np.linspace(0, width/2)
x2 = np.linspace(width/2, width)

plt.rc('text', usetex=True)
plt.rcParams.update({'font.size': 18})
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 4.2))

line1 = ax.plot(x1, psiplus1(x1, 0.3, source, sigmaa), 'g')
line2 = ax.plot(x1, psimin(x1, 0.3, source, sigmaa, width), 'orange')
line3 = ax.plot(x1, psiplus1(x1, 0.7, source, sigmaa), 'r')
line4 = ax.plot(x1, psimin(x1, 0.7, source, sigmaa, width), 'purple')
line5 = ax.plot(x2, psiplus2(x2, 0.3, source, sigmaa, width), 'g')
line6 = ax.plot(x2, psiplus2(x2, 0.7, source, sigmaa, width), 'r')

ax.set_xlim(0-0.01, 1+0.01)
ax.set_xticks([0, width/2, width])
ax.set_xticklabels(['0', 'W/2', 'W'])
ax.legend(('$\psi_+,\ \mu=0.3$', '$\psi_-,\ \mu=0.3$', '$\psi_+,\ \mu=0.7$', '$\psi_-,\ \mu=0.7$'))
plt.xlabel('Depth in Slab (cm)')
plt.ylabel('Angular Flux $\psi($cm^{-2}\ s^{-1}\ ev^{-1}\ rad^{-1}$)$')
plt.tight_layout()
plt.grid()
plt.show()


def phi1(x, source, sigmaa, width):
    return source/(2*sigmaa) * (2-expn(2, sigmaa*(width/2-x))-expn(2, sigmaa*x))


def phi2(x, source, sigmaa, width):
    return source/(2*sigmaa) * (expn(2, sigmaa*(width/2+x)) - expn(2, sigmaa*x))

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 4.2))
ax.plot(x1, phi1(x1, source, sigmaa, width))
ax.plot(x2, phi2(x2, source, sigmaa, width))
plt.grid()
plt.tight_layout()
plt.xlabel('Depth in Slab (cm)')
plt.ylabel('Scalar Flux $\psi($cm^{-2}\ s^{-1}\ ev^{-1}$')
ax.set_xlim(0-0.01, 1+0.01)
ax.set_xticks([0, width/2, width])
ax.set_xticklabels(['0', 'W/2', 'W'])
plt.show()







