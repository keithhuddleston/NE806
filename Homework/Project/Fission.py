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
from Project_Utilities import Nuclide_Data # File written for class
from Removal_Matrix import Interpolate_Group_Total
from Removal_Matrix import Background_Cross_Section
from Scatter_Matrix import Scatter_Matrix
from Watt_Spectrum import Chi_Matrix

# ============================================================================
# Functions
# ============================================================================
def Interpolate_Group_Fission(N, Nuclides, Temperature, Dilution, group_structure):
    # NOTE, I have not included temperature interpolation yet but this will
    # assert statement needs to be here for when this is implemented
    assert 300 <= Temperature <= 1200, 'Temperature value outside of bounds' 
    
    # Because interpolation is not implemented we need to asser this
    assert Temperature in [300, 600, 900, 1200], 'Interpolation not implemented'
    
    # Make sure Dilution is between interprable bounds
    assert 1 <= Dilution <= 10000, 'Background/dilution cross-section out of range'
    
    S = np.zeros((len(group_structure)-1, len(group_structure)-1))
    
    # NOTE this can be buggy later is the calculated dilution cross sections 
    # were not at [1, 10, 100, 1000, 10000] or the values don't follow the 
    # same group structure
    Dilution_Bounds = [1, 10, 100, 1000, 10000]
    for i in range(len(Dilution_Bounds)):
        if Dilution >= Dilution_Bounds[i]:
            index = i
            break
    for i in range(len(Nuclides)):
        if Nuclides[i].B[2] == 1:
            file_name = 'Data/Group/'+Nuclides[i].N+'G_Fission_'+str(Temperature)+'.txt'
            M = np.loadtxt(file_name, skiprows=1, unpack=True)
            row = M[index+1] - (M[index+1]-M[index])\
                     /(Dilution_Bounds[index+1]-Dilution_Bounds[index])\
                     *(Dilution_Bounds[index+1]-Dilution)
            for j in range(len(group_structure)-1):
                S[j] = S[j] + row*N[index]
    return S

def Fission_Matrix(Nuclides, N, Temperature, Dilution, group_structure, nu):
    shape = len(group_structure) - 1
    F = np.zeros((shape, shape))
    F += Interpolate_Group_Fission(N, Nuclides, Temperature, Dilution, group_structure)
    Chi = Chi_Matrix(group_structure/1e6)
    F = F*Chi*nu
    return F

# ============================================================================
# Load Interpolated Interpreted Plotted Data Files and Doppler-Broadened Data
# ============================================================================
H1 = Nuclide_Data('H1', 1.008, [1, 1, 0], 1)
O16 = Nuclide_Data('O16', 15.995, [1, 1, 0], 16)
U235 = Nuclide_Data('U235', 235.044, [1, 1, 1], 235)
U238 = Nuclide_Data('U238', 238.051, [1, 1, 1], 238)

H1.Load_Doppler_Data([600, 900, 1200])
O16.Load_Doppler_Data([600, 900, 1200])
U235.Load_Doppler_Data([600, 900, 1200])
U238.Load_Doppler_Data([600, 900, 1200])

Casmo_16 = np.array([1.00e1,   8.21e-1, 5.53e-3, 4.00e-6, 1.30e-6, 1.15e-6, 
                     1.097e-6, 1.02e-6, 9.71e-7, 8.50e-7, 6.25e-7, 3.50e-7, 
                     2.80e-7, 1.40e-7,  5.80e-8, 3.00e-8, 1.00e-11])*1e6 # eV
Casmo_16 = Casmo_16[::-1]

# Microscopic Dilution/Background Cross Section
Dilution = [1e1, 1e2, 1e3, 1e4, 1e5] # Barns
# ============================================================================
# Testing
# ============================================================================
if __name__ == '__main__':
    nu = 2.54
    # Function Input
    N = [6, 7, 0.1, 1.9]
    s = [20, 10, 10, 0]
    Names = [H1, O16, U235, U238]
    Temperature = 300
    so = Background_Cross_Section(N, s)
    print("For our approximation the dilution cross-section is "+str(so)+'\n')
    
    F = Interpolate_Group_Fission(N, Names, Temperature, so, Casmo_16)
    C = Chi_Matrix(Casmo_16/1e6)
    
    test = Fission_Matrix(Names, N, 300, so, Casmo_16, nu)