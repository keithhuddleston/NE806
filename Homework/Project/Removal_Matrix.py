""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
     
    Create data files for the flux averaged cross sections.
"""

# ============================================================================
# Import statements
# ============================================================================
import numpy as np
import matplotlib.pyplot as plt

# Files written for this  project
from Utilities.Utilities import Nuclide
from Utilities.Utilities import Background_Cross_Section


# ============================================================================
# Functions
# ============================================================================
def Interpolate_Group(Nuclide, Temperature, Dilution, Group_Structure, Type):
    """ Interpolate from the data files in the group folder. """
    print('Interpolating data for '+Nuclide.N+' at '+str(Dilution)+'...\n')
    shape = len(Group_Structure) - 1
    Matrix = np.zeros((shape, shape))

    if Type == 'Fission':
        if Nuclide.B[2] == 0: # If true then there is no fission data.
            return Matrix

    # Data files do not extend past these bounds, future work may extend this.
    assert 300 <= Temperature <= 1200 and 1 <= Dilution <= 10000, \
    'Temperature or Dilution cross-section out of calculated bounds of data'
    
    # Temperature interpolation not yet implemented, hence restricted.
    assert Temperature in [300, 600, 900, 1200], \
    'Temperature interpolation not implemented.'

    # Cross-sections evaluated for these Diluation values.
    Dilutions = [1, 10, 100, 1000, 10000]
    for i in range(len(Dilutions)):
        if Dilution >= Dilutions[i]:
            index = i
    
    file_name = 'Data/Group/'+str(shape)+'/'+Nuclide.N+'G_'+Type+'_'+\
                str(Temperature)+'.txt'

    data = np.loadtxt(file_name, skiprows=1, unpack=True)
    # Linear Interpolation
    row = data[index] + (data[index+1]-data[index]) \
                      / (Dilutions[index+1]-Dilutions[index]) \
                      * (Dilution-Dilutions[index])
    assert Type == 'Fission' or Type == 'Total', 'Wrong File Type'
    if Type == 'Fission':
        for j in range(shape):
            Matrix[j] = Matrix[j] + row
    elif Type == 'Total':
        Matrix = np.diag(row)
        
    return Matrix

def Removal_Matrix(N, Nuclides, Temperature, Dilution, Group_Structure):
    print('='*79+'\n'+'Calculating the removal matrix\n'+'='*79+'\n')
    shape = len(Group_Structure) - 1
    Matrix = np.zeros((shape, shape))
    for i in range(len(Nuclides)):
        Matrix = Matrix + Interpolate_Group(Nuclides[i], Temperature, 
                                            Dilution, Group_Structure, 
                                            'Total') * N[i]

    print('Removal Matrix all done :)...\n')
    
    return Matrix

# ============================================================================
# Testing
# ============================================================================
if __name__ == '__main__':
    print('Testing the file Removal_Matrix.py\n')
    # Function Input
    N = [6, 7, 0.1, 1.9]
    s = [20, 10, 10, 0]
    
    H1 = Nuclide('H1', 1, 1.008, [1, 1, 0])
    H1.Load_Doppler_Data([600, 900, 1200])
    
    O16 = Nuclide('O16', 16, 15.995, [1, 1, 0])
    O16.Load_Doppler_Data([600, 900, 1200])
    
    U235 = Nuclide('U235', 235, 235.044, [1, 1, 1])
    U235.Load_Doppler_Data([600, 900, 1200])
    
    U238 = Nuclide('U238', 238, 238.051, [1, 1, 1])
    U238.Load_Doppler_Data([600, 900, 1200])

    Casmo_16 = np.array([1.00e1,   8.21e-1,  5.53e-3, 4.00e-6, 1.30e-6, 
                         1.15e-6,  1.097e-6, 1.02e-6, 9.71e-7, 8.50e-7, 
                         6.25e-7,  3.50e-7,  2.80e-7, 1.40e-7, 5.80e-8, 
                         3.00e-8,  1.00e-11])*1e6 # eV
    
    Casmo_2 = np.array([1.00e1, 1.00e-7, 1.00e-11])*1e6 # eV
    
    Nuclides = [H1, O16, U235, U238]
    
    Temperature = 300
    
    so = Background_Cross_Section(N, s)
    
    print("For our approximation the dilution cross-section is "+str(so)+'\n')
    
    Test_1 = Interpolate_Group(H1, Temperature, so, Casmo_16, 'Total')
    
    Test_2 = Removal_Matrix(N, Nuclides, 300, so, Casmo_16)
    
    # TODO fix this error, is messed up because of currenlty used naming conv.
    # Test_4 = Removal_Matrix(N, Nuclides, 300, so, Casmo_2)