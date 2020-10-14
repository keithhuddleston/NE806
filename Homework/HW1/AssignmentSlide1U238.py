# Module Imports =============================================================
import numpy as np
import matplotlib.pyplot as plt

# Define Parameters Used in Breit-Wigner Formula =============================
sigma_p = 11.2934  # potential cross section (in barns)
I = 0  # s wave
E0, J, GN, GG, GFA, GFB = np.loadtxt('u238.txt', unpack=True)
A = 238.0
E = np.logspace(-1, 3, 10000)
    
def breitWigner(E0, J, GN, GG, GFA, GFB, A, E):
    G = GN+GG
    sigma_g = np.zeros_like(E)
    sigma_n = np.zeros_like(E)
    sigma_n[:] = sigma_p
    sigma_e = np.zeros_like(E)
    #  Nuclear Radius in (cm), units cancel out with scattering interfence term
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
    return sigma_e, sigma_g

sigma_e, sigma_g = breitWigner(E0, J, GN, GG, GFA, GFB, A, E)
# E2 = tabulated energy values coresponding to exs, elastic cross section
E2, exs = np.loadtxt('u238e.txt', skiprows=1, delimiter=',', unpack=True)
E3, rcxs = np.loadtxt('u238rc.txt', skiprows=1, delimiter=',', unpack=True)

    
# Plotting ===================================================================
plt.loglog(E, sigma_e)
plt.loglog(E2, exs), plt.xlim(0.1, 1000)
plt.legend(['sigma_e', 'BNL elastic xs']), plt.tight_layout()
plt.xlim(5, 10)
plt.show()

plt.loglog(E, sigma_g)
plt.loglog(E3, rcxs), plt.xlim(0.1, 1000)
plt.legend(['sigma_g', 'BNL radiative capture xs']), plt.tight_layout()
plt.xlim(5, 10)
plt.show()