""" Neutronics Pre-Lecture to-do
"""

#  Imports ===================================================================
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

#  Bring in cross section data ===============================================
E1, U238_ng = np.loadtxt('Data/U238_NG.txt', skiprows=1, unpack=True, delimiter=',')
E2, U238_es = np.loadtxt('Data/U238_NG.txt', skiprows=17, unpack=True, delimiter=',')
E3, H1_es = np.loadtxt('Data/H1_ES.txt', skiprows=1, unpack=True, delimiter=',')

Emesh = np.linspace(0.1, 50, 1000)
U238_ng = np.interp(Emesh, E1, U238_ng)
U238_es = np.interp(Emesh, E2, U238_es)
H1_es = np.interp(Emesh, E3, H1_es)
U238_tot = U238_ng + U238_es
plt.rcParams.update({'font.size': 18})
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
ax.loglog(Emesh, U238_ng)
ax.loglog(Emesh, U238_es)
ax.loglog(Emesh, H1_es)
plt.tight_layout(), plt.grid()
plt.xlabel('Energy [eV]'), plt.ylabel('$\sigma(E) [Barns]$')
plt.show()

# ============================================================================
def alpha(A):
    return ((A-1)/(A+1))**2

def scatter_probability(E_i, E_j, alpha) :
    p = (1.0/E_j/(1.0-alpha)) * 1.0*((E_i >= alpha*E_j))
    return p 


def compute_spectrum(E, Sigma_t, Sigma_s, N, alpha) :
    N_E = len(E)
    phi = np.zeros(N_E)
    phi[N_E-1] = 1.0
    for i in range(N_E-2, -1, -1) :
        Q_i = 0.0
        for j in range(N_E-1, i, -1) :
            dE = E[j] - E[j-1]
            E_bar = np.sqrt(E[j]*E[j-1])
            Q_i += phi[j] * (Sigma_s[j]*N) * scatter_probability(E[i], E[j], alpha) * dE
        phi[i] = Q_i / Sigma_t[i]
    return phi

# Functions ==================================================================
def phinr(sd, st, E):
    return 1 / (E*(st+sd))

def phiwr(sd, sa, E):
    return 1 / (E*(sa+sd))

sd = 100*H1_es
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
ax.loglog(Emesh, phinr(sd,U238_tot,Emesh)/phinr(sd,U238_tot,Emesh)[-1])
ax.loglog(Emesh, phiwr(sd,U238_ng,Emesh)/phiwr(sd,U238_ng,Emesh)[-1])
ax.loglog(Emesh, 1/Emesh*Emesh[-1])
plt.xlabel('Energy [eV]'), plt.ylabel('Cross Section [Barns]')
plt.legend(['$\phi_{nr}$','$\phi_{wr}$','1/E'])
plt.tight_layout(), plt.grid()

sigt = 1.0*(U238_tot)+100.0*H1_es
phi1 = compute_spectrum(Emesh, sigt, U238_es, 1, alpha(238)) 
phi2 = compute_spectrum(Emesh, sigt, H1_es, 100, alpha(1)) 
phi = phi1+phi2
ax.loglog(Emesh, phi)





















