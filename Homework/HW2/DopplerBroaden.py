""" Neutronics HW2
"""
# Module Imports =============================================================
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.special import erf
sp.init_printing(use_unicode=True)
plt.rcParams.update({'font.size': 18})

# Define Parameters Used in Breit-Wigner Formula =============================
sigma_p = 11.2934  # potential cross section (in barns)
I = 0  # s wave
E0, J, GN, GG, _, _ = np.loadtxt('u238.txt', unpack=True)
G = GN+GG
A = 238.0
E = np.logspace(-1, 3, 10000)
sigma_g = np.zeros_like(E)
sigma_n = np.zeros_like(E)
sigma_n[:] = sigma_p

sigma_e = np.zeros_like(E)

#  Nuclear Radius in (cm)
R = 1.25*10**-13 * 238**(1/3)

for i in range(len(E0)):
    # statistical spin factor, where I=nuclear spint, J=total spin
    g = (2*J[i]+1)/(2*(2*I+1))

    # total cross section at resonance energy (DH 2-36)
    sigma_0 = 2.608e6 * (A+1)**2/(A**2 * E0[i]) * (GN[i]/G[i]) * g

    # capture cross section (DH 2-35)
    y = (2/G[i])*(E - E0[i])
    sigma_g += sigma_0 * (GG[i]/G[i]) * np.sqrt(E0[i]/E) * (1/(1+y**2))
    sigma_e += sigma_0 * (GN[i]/G[i]) * np.sqrt(E0[i]/E) * (1/(1+y**2)) \
             + sigma_0 * (2*R/(4.55*10**-10/E0[i]**0.5))*(y/(1+y**2))

sigma_e = sigma_e + sigma_n  # Add potential scattering

# Plot Breit-Wivner ==========================================================
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
ax.loglog(E, sigma_g, 'orange')
ax.set_xlabel('Neutron Energy (eV)'), ax.set_ylabel('Cross Section (Barns)')
ax.legend(['$Breit-Wigner\ \sigma_\gamma$', '$BNL\ \sigma_\gamma$'], loc=3),\
          plt.tight_layout()
ax.grid(), plt.show()

#  Imports ===================================================================
F0 = lambda a: erf(a)
H0 = lambda a, b: F0(a) - F0(b)

F1 = lambda a: np.sqrt(1/np.pi) * (1-np.exp(-a**2))
H1 = lambda a, b: F1(a) - F1(b)

F2 = lambda a: (1/2)*erf(a) - (a/np.sqrt(np.pi))*np.exp(-a**2)
H2 = lambda a, b: F2(a) - F2(b)

F3 = lambda a: np.sqrt(1/np.pi) * (1-(1+a**2)*np.exp(-a**2))
H3 = lambda a, b: F3(a) - F3(b)

F4 = lambda a: (3/4)*erf(a) - np.sqrt(1/np.pi)*((3*a/2)+a**3)*np.exp(-a**2)
H4 = lambda a, b: F4(a) - F4(b)


def Akf(E1, E2, S1, S2):
    den = (E2 - E1)
    num = (E2*S1) - (E1*S2)
    return num/den


def Ckf(E1, E2, S1, S2, alpha):
    den = (E2 - E1)*alpha
    num = (S2 - S1)
    return num/den


def doppler(E2, E1, S1, T1, T2, M, m=1.009):
    """

    Parameters
    ----------
    E2 : array_like
        New energy mesh at which to evaluate cross sections at.
    E1 : array_like
        Energy mesh of the reference cross section data.
    S1 : array_like
        Cross section data of reference.
    T1 : Float
        Temperature of the reference cross section data.
    T2 : Float
        New temperature at which to evaluate cross sections at.
    m : Float
        Atomic mass of projectile, 1.009 for Neutron.
    M : Float
        Atomic mass of target.
    neg : Int, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    S2 : array_like
        Reavaluated cross section data for energies E2, and temperature T2

    """
    Bk     = 6.617*10**-5        # Boltzman Constant, [ev K^-1]
    alpha  = (M/m)/(Bk*(T2-T1))  # Alpha term found in Doppler broadening Eqs.
    S2     = np.zeros(len(E2))   # Initialize new cross section data array
    S2_pos = np.zeros(len(E2))
    S2_neg = np.zeros(len(E2))
    
    # Evaluate Doppler-broadened cross section at specified energy E2[i]
    for i in range(len(E2)):
        S2i = 0
        y = [-1*np.sqrt(alpha*E2[i]), np.sqrt(alpha*E2[i])]
        for j in range(len(y)):
            Ek1 = E1[:-1]
            Ek2 = E1[1:]
            Sk1 = S1[:-1]
            Sk2 = S1[1:]            
            xk1 = np.sqrt(alpha*Ek1)
            xk2 = np.sqrt(alpha*Ek2)

            Ak = Akf(Ek1, Ek2, Sk1, Sk2)
            Ck = Ckf(Ek1, Ek2, Sk1, Sk2, alpha)

            Zk1 = xk1 - y[j]
            Zk2 = xk2 - y[j]

            H0k = H0(Zk2, Zk1)
            H1k = H1(Zk2, Zk1)
            H2k = H2(Zk2, Zk1)
            H3k = H3(Zk2, Zk1)
            H4k = H4(Zk2, Zk1)

            S2i = H4k * (Ck) \
                 + H3k * (4*Ck*y[j]) \
                 + H2k * (Ak+6*Ck*y[j]**2) \
                 + H1k * (2*Ak*y[j]+4*Ck*y[j]**3) \
                 + H0k * (Ak*y[j]**2+Ck*y[j]**4)
            S2i = sum(S2i)
            if j == 0:
                S2_neg[i] = S2i/2/y[j]**2
            else:
                S2_pos[i] = S2i/2/y[j]**2
            S2 = S2_pos - S2_neg
    return S2

M = 238.051
Emesh = np.logspace(-5, 7, 10000)
y1 = doppler(Emesh, E, sigma_g, 0, 300, M)
y2 = doppler(Emesh, E, sigma_g, 0, 600, M)
y3 = doppler(Emesh, E, sigma_g, 0, 900, M)
np.savetxt('U238_ES_0_BW.txt', np.vstack(np.transpose([E, sigma_g])), delimiter=',')
np.savetxt('U238_ES_300_BW.txt', np.vstack(np.transpose([Emesh, y1])), delimiter=',')
np.savetxt('U238_ES_600_BW.txt', np.vstack(np.transpose([Emesh, y2])), delimiter=',')
np.savetxt('U238_ES_900_BW.txt', np.vstack(np.transpose([Emesh, y3])), delimiter=',')
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
ax.loglog(E, sigma_g, Emesh, y1, Emesh, y2, Emesh, y3), ax.set_xlim(1, 1000)
plt.legend(['$Breit-Wigner\ \sigma_\gamma, T_2=0K$', '$T_2 = 300K$', '$T_2 = 600K$', '$T_2 = 900K$'])



