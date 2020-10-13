""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
"""

# ============================================================================
# Import statements
# ============================================================================
import numpy as np
import matplotlib.pyplot as plt
from NE806_Functions import doppler # File written for class

# ============================================================================
# Load in BNL interpreted data files from the Data subfolder
# ============================================================================
# Casmo 2 group structure
c2 = np.array([.001e1, 1.00e-7, 1.00e-11])*1e6 # eV

# Casmo 16 group structure
c16 = np.array([1.00e1 , 8.21e-1, 5.53e-3, 4.00e-6, 1.30e-6, 1.15e-6, 1.097e-6,
                1.02e-6, 9.71e-7, 8.50e-7, 6.25e-7, 3.50e-7, 2.80e-7, 1.40e-7 ,
                5.80e-8, 3.00e-8, 1.00e-11])*1e6 # eV

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
    def __init__(self, M):
        self.M = M # Atomic Mass
        
        self.file_contents = ['ES', 'NG']
        self.base_temp = 0 # K
        
        self.ESXS = {}
        self.ESEM = {}
        self.NGXS = {}
        self.NGEM = {}
        return
    
    def load_plot_data(self, filename):
        """ Loads BNL interpreted plotting data
        """
        ES, XS = np.loadtxt(filename, delimiter=',', unpack=True, skiprows=1)
        return ES, XS
    
    def load_data(self, filename, file_content, temperature):
        EM, XS = np.loadtxt(filename, delimiter=',', unpack=True, skiprows=1)
        if file_content == 'ES':
            self.ESXS[str(temperature)] = XS
            self.ESEM[str(temperature)] = EM
        elif file_content == 'NG':
            self.NGXS[str(temperature)] = XS
            self.NGEM[str(temperature)] = EM
      
# ============================================================================
# Testing Oct 12. 2020
# ============================================================================

# Load and Doppler Broaden data here
Flag = True
if Flag:
    Temps = [600, 900, 1200] # Temperature values specified in HW
    
    H1 = isotope(M=1.008)
    H1.load_data('Data/H1_ES.txt', 'ES', 300)
    
    O16 = isotope(M=15.995)
    O16.load_data('Data/O16_ES.txt', 'ES', 300)
    
    U235 = isotope(M=235.044)
    U235.load_data('Data/U235_ES.txt', 'ES', 300)
    U235.load_data('Data/U235_NG.txt', 'NG', 300)
    
    U238 = isotope(M=238.051)
    U238.load_data('Data/U238_ES.txt', 'ES', 300)
    U238.load_data('Data/U238_NG.txt', 'NG', 300)
    
    # Emesh based on Casmo group structure
    Emesh = np.logspace(-5, 7, 75000)
    
    # Doppler Broaden H1
    for i in Temps:
        E1 = H1.ESEM['300']
        XS = H1.ESXS['300']
        H1.ESXS[str(i)] = doppler(Emesh, E1, XS, 300, i, H1.M)
        H1.ESEM[str(i)] = Emesh
        
    # Doppler Broaden O16
    for i in Temps:
        E1 = O16.ESEM['300']
        XS = O16.ESXS['300']
        O16.ESXS[str(i)] = doppler(Emesh, E1, XS, 300, i, H1.M)
        O16.ESEM[str(i)] = Emesh

plt.rcParams.update({'font.size': 18})
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
for key in H1.ESEM.keys():
    ax.loglog(H1.ESEM[key], H1.ESXS[key])
plt.tight_layout(), plt.legend(H1.ESXS.keys())
plt.show()

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
for key in O16.ESEM.keys():
    ax.loglog(O16.ESEM[key], O16.ESXS[key])
plt.tight_layout(), plt.legend(O16.ESXS.keys())
plt.show()























