""" NE 806, Neutronics, Function to recontruct cross section data using the 
    Breit Wigner method

    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
"""

# Module Imports =============================================================
import numpy as np
import matplotlib.pyplot as plt

# Define Parameters Used in Breit-Wigner Formula =============================
sigma_p = 11.2934  # potential cross section (in barns)
A = 238.0  # Mass Number
I = 0  # s wave

# Load in resonance parameters
E0, J, GN, GG, GFA, GFB = np.loadtxt('u238.txt', unpack=True)
G = GN + GG
# Define Energy Mesh to reconstruct cross section data
E = np.linspace(0.1, 1000, 50000)

sigma_g = np.zeros_like(E)  # Neutron Gamma interaction (Radiative Capture) 
sigma_e = np.zeros_like(E)  # Elastic Scatter

#  Nuclear Radius in (cm)
R = 1.25*10**-13 * 238**(1/3) # Nuclear radius approximation

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

sigma_e = sigma_e + sigma_p  # Add potential scattering

# Plotting ===================================================================
fig, ((ax1), (ax2)) = plt.subplots(nrows=2, ncols=1, figsize=(12,12))
ax1.loglog(E, sigma_g, 'orange')
ax1.set_xlabel('Neutron Energy (eV)')
ax1.set_ylabel('Cross Section (Barns)')
ax1.legend(['$Breit-Wigner\ \sigma_\gamma$'], loc=3), plt.tight_layout()
ax2.loglog(E, sigma_e, 'red')
ax2.set_xlim(0.1, 1000)
ax2.set_xlabel('Neutron Energy (eV)')
ax2.set_ylabel('Cross Section (Barns)')
ax2.legend(['$Breit-Wigner\ \sigma_e$'], loc=3), plt.tight_layout()
plt.show()