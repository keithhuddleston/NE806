""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
"""

# ============================================================================
# Import statements
# ============================================================================
import matplotlib.pyplot as plt
import numpy as np


# ============================================================================
# Functions
# ============================================================================
def Watt_Spectrum(E):
    return 0.4865*np.sinh(np.sqrt(2*E))*np.exp(-E)

def Chi_Matrix(group_structure):
    shape = len(group_structure) - 1
    energy_group = []
    chi_group = []
    for i in range(shape):
        start = group_structure[i]
        end = group_structure[i+1]
        energy = np.linspace(start, end, 10)
        energy_group.append(energy)
        chi_group.append(Watt_Spectrum(energy))
    for i in range(shape):
        chi_group[i] = np.trapz(chi_group[i], energy_group[i])
    chi_matrix = np.zeros((shape, shape))
    for i in range(shape):
        chi_matrix[i] = np.zeros(shape)+chi_group[i]
    return chi_matrix


# ============================================================================
# Testing
# ============================================================================
Casmo_16 = np.array([1.00e1,   8.21e-1, 5.53e-3, 4.00e-6, 1.30e-6, 1.15e-6, 
                     1.097e-6, 1.02e-6, 9.71e-7, 8.50e-7, 6.25e-7, 3.50e-7, 
                     2.80e-7, 1.40e-7,  5.80e-8, 3.00e-8, 1.00e-11]) # MeV
Casmo_16 = Casmo_16[::-1]

# Microscopic Dilution/Background Cross Section
Dilution = [1e1, 1e2, 1e3, 1e4, 1e5] # Barns

if __name__ == '__main__':
    Chi = Chi_Matrix(Casmo_16)
    
    # TODO add some dope plotting