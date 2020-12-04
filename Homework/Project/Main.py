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
from Project_Utilities import Nuclide_Data
from Removal_Matrix import Background_Cross_Section
from Removal_Matrix import Removal_Matrix
from Scatter_Matrix import Scatter_Matrix
from Fission_Matrix import Fission_Matrix
from scipy.linalg import eig

# ============================================================================
# Load Interpolated Interpreted Plotted Data Files and Doppler-Broadened Data
# ============================================================================
H1 = Nuclide_Data('H1', 1.008, [1, 1, 0], 1)
H1.Load_Doppler_Data([600, 900, 1200])

O16 = Nuclide_Data('O16', 15.995, [1, 1, 0], 16)
O16.Load_Doppler_Data([600, 900, 1200])

U235 = Nuclide_Data('U235', 235.044, [1, 1, 1], 235)
U235.Load_Doppler_Data([600, 900, 1200])

U238 = Nuclide_Data('U238', 238.051, [1, 1, 1], 238)
U238.Load_Doppler_Data([600, 900, 1200])

Casmo_16 = np.array([1.00e1,   8.21e-1, 5.53e-3, 4.00e-6, 1.30e-6, 1.15e-6, 
                     1.097e-6, 1.02e-6, 9.71e-7, 8.50e-7, 6.25e-7, 3.50e-7, 
                     2.80e-7,  1.40e-7, 5.80e-8, 3.00e-8, 1.00e-11])*1e6 # eV

Casmo_16 = Casmo_16[::-1]

# Microscopic Dilution/Background Cross Section
Dilution = [1e1, 1e2, 1e3, 1e4, 1e5] # Barns
# ============================================================================
# Testing
# ============================================================================
if __name__ == '__main__':
    N = [6, 7, 0.1, 1.9]
    s = [20, 10, 10, 0]
    Names = ['H1', 'O16', 'U235', 'U238']
    Temperature = 300
    so = Background_Cross_Section(N, s)
    print("For our approximation the dilution cross-section is "+str(so)+'\n')  
    T = Interpolate_Group_Total(N, Names, Temperature, so, Casmo_16)
    Nuclides = [H1, O16, U235, U238]
    print('Calculating S this takes awhile...\n')
    S = Scatter_Matrix(Nuclides, N, Casmo_16, so, 300)
    
    R = np.T - S
    Names = [H1, O16, U235, U238]
    nu = 2.54
    F = Fission_Matrix(Names, N, 300, so, Casmo_16, nu)
    
    test = eig(R, F)