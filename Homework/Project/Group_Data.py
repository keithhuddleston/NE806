""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
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

def Narrow_Resonance_Approx(xs_total, xs_dilution, energy_total):
    return 1/(energy_total*(xs_total + xs_dilution))

# Microscopic Flux Weighted Average Cross Section.  Note, that the energy mesh
# must be ultra fine for this simplifaction to work properly.
def sigma_g(Phi, XS):
    return(sum(Phi*XS)/sum(Phi))


def Seperate(x_vals, y_vals, Group_Structure):
    """ Seperate input data into groups based on given group structure """
    
    index = []
    step = 0
    for x in x_vals:
        if x < Group_Structure[step]:
            step = step + 1
            index.append
    return index

def Make_Group_Data(N, Group_Structure, Dilution):
    """ Create flux weighted averages cross section data """
    Temperatures = N.T
    for t in Temperatures:
        for i in N.B:
            if i:
                # Index '1' should always correspond to total cross section
                phi = Narrow_Resonance_Approx(N.XS[t][1], Dilution[0], 
                                              N.EM[t][1], N.EM[t][i])
                phi = Seperate(N.EM[t][i], phi, Group_Structure)
                xs = Seperate(N.EM[t][i], N.XS[t][i])
            xs_group = np.zeros(len(Group_Structure))
            for j in range(len(Group_Structure)):
                xs_group[j] = sigma_g(phi[j], xs[j])
    return


# ============================================================================
# Define Group and Background Data
# ============================================================================
Casmo_2 = np.array([1.00e1, 1.00e-7, 1.00e-11])*1e6 # eV
Casmo_2 = Casmo_2[::-1]

Casmo_16 = np.array([1.00e1,   8.21e-1, 5.53e-3, 4.00e-6, 1.30e-6, 1.15e-6, 
                     1.097e-6, 1.02e-6, 9.71e-7, 8.50e-7, 6.25e-7, 3.50e-7, 
                     2.80e-7, 1.40e-7,  5.80e-8, 3.00e-8, 1.00e-11])*1e6 # eV
Casmo_16 = Casmo_16[::-1]

# Microscopic Dilution/Background Cross Section
Dilution = [1e1, 1e2, 1e3, 1e4, 1e5]


# ============================================================================
# Load Interpolated Interpreted Plotted Data Files and Doppler-Broadened Data
# ============================================================================
H1 = Nuclide_Data('H1', 1.008, [1, 1, 0])
# O16 = Nuclide_Data('O16', 15.995, [1, 1, 0])
# U235 = Nuclide_Data('U235', 235.044, [1, 1, 1])
# U238 = Nuclide_Data('U238', 238.051, [1, 1, 1])

H1.Load_Doppler_Data([600, 900, 1200])
# O16.Load_Doppler_Data([600, 900, 1200])
# U235.Load_Doppler_Data([600, 900, 1200])
# U238.Load_Doppler_Data([600, 900, 1200])


# ============================================================================
# Create Group Files
# ============================================================================
Make_Group_Data(H1, Casmo_16, Dilution)
