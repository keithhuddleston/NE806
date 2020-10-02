# Module Imports =============================================================
import numpy as np
import matplotlib.pyplot as plt

# Define Parameters Used in Breit-Wigner Formula =============================
sigma_p = 11.5860  # potential cross section (in barns)
I = 0  # s wave
E0, J, GN, GG, GFA, GFB = np.loadtxt('u235processed.txt', unpack=True)
G = GN+GG
A = 235.0
E = np.logspace(-1, 3, 10000)
sigma_f = np.zeros_like(E)

for i in range(len(E0)):
    # statistical spin factor, where I=nuclear spint, J=total spin
    g = (2*J[i]+1)/(2*(2*I+1))

    # total cross section at resonance energy (DH 2-36)
    sigma_0 = 2.608e6 * (A+1)**2/(A**2 * E0[i]) * (GN[i]/G[i]) * g

    # capture cross section (DH 2-35)
    y = (2/G[i])*(E - E0[i])
    sigma_f += sigma_0 * (GFA[i]/G[i]) * np.sqrt(E0[i]/E) * (1/(1+y**2))


# E2 = tabulated energy values coresponding to exs, elastic cross section
E2, exs = np.loadtxt('u238e.txt', skiprows=1, delimiter=',', unpack=True)

# Plotting ===================================================================
plt.loglog(E, sigma_f)
plt.show()