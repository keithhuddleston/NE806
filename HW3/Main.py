""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
"""

# ============================================================================
# Import statements
# ============================================================================
import numpy as np
import matplotlib.pyplot as plt
import re

# ============================================================================
# Load in BNL interpreted data files from the Data subfolder
# ============================================================================
# Casmo 2 group structure
c2 = np.array([.001e1, 1.00e-7, 1.00e-11])*1e6

# Casmo 16 group structure
c16 = np.array([1.00e1 , 8.21e-1, 5.53e-3, 4.00e-6, 1.30e-6, 1.15e-6, 1.097e-6,
       1.02e-6, 9.71e-7, 8.50e-7, 6.25e-7, 3.50e-7, 2.80e-7, 1.40e-7 ,
       5.80e-8, 3.00e-8, 1.00e-11])*1e6

def load_BNL_RP(filename):
    """ File for loading BNL interpreted resonance parameter files
    """
    file = open(filename, 'r')
    lines = file.readlines()
    l = len(lines)
    data = np.zeros((l, 6)) # There should always be six columns in this file
    
    # Alter file values so that python can recognize them as floating point
    p1 = r'(\+*\-*\d\.\d*\+*\-*\d+)'
    for i in range(l):
        line = lines[i]
        found = re.findall(p1, line)
        for j in range(6):
            value = found[j]
            for k in range(len(value))[::-1]:
                if value[k] == '-' or value[k] == '+':
                    value = value[:k]+'e'+value[k:]
                    break
            data[i][j] = value
            
    # Seperate data into column vectors containing E0, J, GN, GG, GFA, GFB
    E0, J, GN, GG, GFA, GFB = [np.zeros(l) for i in range(6)]
    for i in range(l):
        E0[i], J[i], GN[i], GG[i], GFA[i], GFB[i] = data[i]
        ind = 0
    for i in range(len(E0)):
        if E0[i] < 0:
            ind = i+1
    return E0[ind:], J[ind:], GN[ind:], GG[ind:], GFA[ind:], GFB[ind:]

def breitWigner(E0, J, GN, GG, GFA, GFB, A, sigma_p, I, E):
    """ Reconstructs cross section data given resonance parameters
    """
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
def fwxs(E, X, Phi, GS):
    """
    Input:
          E, energy mesh of X and Phi
          X, Cross section data
          Phi, flux
          GS, specified group structure
    Output:
          GX,
    """
    # assert len(E) == len(X) == len(Phi), 'Length of input must be the same'
    # Step 1 add e values of GS to E
    # Step 2 interpolate values for new Es
    # Step 3 perform integral
    plt.loglog(E, X)
    plt.vlines(GS, min(X), max(X), linestyle='--'), plt.show()

# ============================================================================
# Classes
# ============================================================================
class isotope:
    def __init__(self):
        return
    
    def addRP(self, E0, J, GN, GG, GFA, GFB, ):
        self.E0 = E0
        self.J = J
        self.GN = GN
        self.GG = GG 
        self.GFA = GFA
        self.GFB = GFB
      
# ============================================================================
# Testing
# ============================================================================

# Test 1. Check that function form of BW is working right
E, XS = np.loadtxt('Data/U238_NG.txt', delimiter=',', unpack=True, skiprows=1)

I = 0 
A = 238
sigma_p = 11.2934
E1 = np.logspace(-11, 10, 50000)
E0, J, GN, GG, GFA, GFB = load_BNL_RP('Data/U238_RP_L0.txt')
sigma_e, sigma_g = breitWigner(E0, J, GN, GG, GFA, GFB, A, sigma_p, I, E1)

# Compare Dr. Robert's resonance file to full file
E01, J1, GN1, GG1, GFA1, GFB1 = np.loadtxt('u238.txt', unpack=True)

# Resonance recontruction using reduced file
sigma_e0, sigma_g0 = breitWigner(E01, J1, GN1, GG1, GFA1, 
                                 GFB1, A, sigma_p, I, E1)

# Do again for l=1 resonance data
E0, J, GN, GG, GFA, GFB = load_BNL_RP('Data/U238_RP_L1.txt')
sigma_e1, sigma_g1 = breitWigner(E0, J, GN, GG, GFA, GFB, A, sigma_p, I, E1)

# Plot results
fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(12,24))
ax1.loglog(E, XS, E1, sigma_g0, E1, sigma_g, E1, sigma_g1)
ax1.legend(['BNL interpreted', 'Reduced l=0', 'Full l=0', 'Full l=1'])
ax1.set_xlim(5, 10)
ax1.set_ylabel('Cross Section (Barns)')

ax2.loglog(E, XS, E1, sigma_g0, E1, sigma_g, E1, sigma_g1)
ax2.set_xlim(700, 1000)
ax2.set_ylabel('Cross Section (Barns)')

ax3.loglog(E, XS, E1, sigma_g0, E1, sigma_g, E1, sigma_g1)
ax3.set_xlim(1000, 1300)
ax3.set_ylabel('Cross Section (Barns)')

ax4.loglog(E, XS, E1, sigma_g0, E1, sigma_g, E1, sigma_g1)
ax4.vlines(c16[1:], ymin=min(sigma_e), ymax=max(sigma_e), linestyle='--')
ax4.set_xlim(min(c16), max(c16))
ax4.set_xlabel('Energy (eV)')
ax4.set_ylabel('Cross Section (Barns)')
plt.show()