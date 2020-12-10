""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
     
    Create data files for the flux averaged cross sections. Note that plotting
    is being performed in this file as well, in contrast to the Doppler-Broad-
    ening step.
"""

# ============================================================================
# Import statements
# ============================================================================
import numpy as np
import matplotlib.pyplot as plt

# File written for class
from Utilities.Utilities import Nuclide
from Utilities.Utilities import Background_Cross_Section
from Removal_Matrix import Removal_Matrix
from Scatter_Matrix import Scatter_Matrix
from Fission_Matrix import Fission_Matrix

# ============================================================================
# Load Interpolated Interpreted Plotted Data Files and Doppler-Broadened Data
# ============================================================================
H1 = Nuclide('H1', 1, 1.008, [1, 1, 0])
H1.Load_Doppler_Data([600, 900, 1200])

O16 = Nuclide('O16', 16, 15.995, [1, 1, 0])
O16.Load_Doppler_Data([600, 900, 1200])

U235 = Nuclide('U235', 235, 235.044, [1, 1, 1])
U235.Load_Doppler_Data([600, 900, 1200])

U238 = Nuclide('U238', 238, 238.051, [1, 1, 1])
U238.Load_Doppler_Data([600, 900, 1200])

# ============================================================================
# Testing
# ============================================================================
if True:
    Nuclides = [H1, O16, U235, U238]
    N = [6, 7, 10, 1.9]
    s = [20, 10, 10, 0]
    Temperature = 300
    nu = 2.54
    K = []
    Test = np.linspace(0.1, 10, 10)
    for i in Test:
        N[2] = i
        so = Background_Cross_Section(N, s)
        print("For our approximation the dilution cross-section is "+str(so)+'\n')  
    
    
        print('Calculating S this takes awhile...\n')
        Casmo_16 = np.array([1.00e1,   8.21e-1,  5.53e-3, 4.00e-6, 1.30e-6, 
                             1.15e-6,  1.097e-6, 1.02e-6, 9.71e-7, 8.50e-7, 
                             6.25e-7,  3.50e-7,  2.80e-7, 1.40e-7, 5.80e-8, 
                             3.00e-8,  1.00e-11])*1e6 # eV
        Casmo_2 = np.array([1.00e1, 1.00e-7, 1.00e-11])*1e6 # eV
        R = Removal_Matrix(N, Nuclides, 300, so, Casmo_2)
        Casmo_16 = Casmo_16[::-1]
        Casmo_2 = Casmo_2[::-1]
        S = Scatter_Matrix(Nuclides, N, Casmo_2, so, 300)
        T = R - S
        plt.matshow(R)
        plt.matshow(S)
        plt.matshow(T)
        Casmo_16 = np.array([1.00e1,   8.21e-1,  5.53e-3, 4.00e-6, 1.30e-6, 
                             1.15e-6,  1.097e-6, 1.02e-6, 9.71e-7, 8.50e-7, 
                             6.25e-7,  3.50e-7,  2.80e-7, 1.40e-7, 5.80e-8, 
                             3.00e-8,  1.00e-11])  # MeV
        Casmo_2 = np.array([1.00e1, 1.00e-7, 1.00e-11])*1e6 # eV
        F = Fission_Matrix(nu, N, Nuclides, 300, so, Casmo_2)
        
        A = np.matmul(np.linalg.inv(T), F)
        k = np.linalg.eig(A)
        K.append(k)
