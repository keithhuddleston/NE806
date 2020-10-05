""" NE 806, Neutronics Homework 3

    Author: Keith Huddleston
    Email 1: kdhuddle@ksu.edu
    Email 2: keithhuddle@gmail.com
"""

# Module Imports =============================================================
import numpy as np
import openmc as mc

# Group Structure ============================================================

# Casmo 2 group structure
c2 = [.001e1, 1.00e-7, 1.00e-11]

# Casmo 16 group structure
c16 = [1.00e1 , 8.21e-1, 5.53e-3, 4.00e-6, 1.30e-6, 1.15e-6, 1.097e-6,
       1.02e-6, 9.71e-7, 8.50e-7, 6.25e-7, 3.50e-7, 2.80e-7, 1.40e-7 ,
       5.80e-8, 3.00e-8, 1.00e-11]

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

def phinr(sd, st, E):
    return 1 / (E*(st+sd))

def phiwr(sd, sa, E):
    return 1 / (E*(sa+sd))

def fwxs(E, sigma, group):
    
    return np.interp(group, E, sigma)