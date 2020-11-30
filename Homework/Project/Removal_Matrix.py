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


# ============================================================================
# Functions
# ============================================================================
def Background_Cross_Section(N, s):
    """ 
    Calculates background cross section for the case that H1, O16, U235, and
    U238 are the only nuclides that the core is composed of.  NOTE this is it-
    self an approximation.
    
    N = [N_H1, N_O16, N_U235, N_U238]
    S = [s_H1_es, s_O16_es, s_U235_potential, s_U238 not needed actually]
    """
    # Hmm... actually we only need three of the cross sections so will leave 
    # the last value in "s" as none for now.
    assert len(N) == len(s), \
    'There must be a number density "N" for each cross-section "s".'
    
    num = N[0]*s[0] + N[1]*s[1] + N[2]*s[2]
    den = N[3]
    return num/den

def Interpolate_Group_Total(N, Names, Temperature, Dilution, group_structure):
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
    for i in range(len(Names)):
        file_name = 'Data/Group/'+Names[i]+'G_Total_'+str(Temperature)+'.txt'
        M = np.loadtxt(file_name, skiprows=1, unpack=True)
        diagonal = M[index+1] - (M[index+1]-M[index])\
                 /(Dilution_Bounds[index+1]-Dilution_Bounds[index])\
                 *(Dilution_Bounds[index+1]-Dilution)
        S = S + np.diag(diagonal)*N[index]
    return S


# ============================================================================
# Testing
# ============================================================================
Casmo_16 = np.array([1.00e1,   8.21e-1, 5.53e-3, 4.00e-6, 1.30e-6, 1.15e-6, 
                     1.097e-6, 1.02e-6, 9.71e-7, 8.50e-7, 6.25e-7, 3.50e-7, 
                     2.80e-7, 1.40e-7,  5.80e-8, 3.00e-8, 1.00e-11])*1e6 # eV
if __name__ == '__main__':
    # Function Input
    N = [6, 7, 0.1, 1.9]
    s = [20, 10, 10, 0]
    Names = ['H1', 'O16', 'U235', 'U238']
    Temperature = 300
    so = Background_Cross_Section(N, s)
    print("For our approximation the dilution cross-section is "+str(so)+'\n')
    
    R = Interpolate_Group_Total(N, Names, Temperature, so, Casmo_16)