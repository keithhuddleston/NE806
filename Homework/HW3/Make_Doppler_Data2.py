""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
     
    Note, expect this file to run for awhile!
"""

# ============================================================================
# Import statements
# ============================================================================
import numpy as np
from NE806_Functions import Doppler # File written for class
from Plot_BNL_Data import Nuclide_Data

# ============================================================================
# Create Text Files Containing the Doppler Broadened Cross Sections
# ============================================================================

# Note, the data we are performing Doppler Broadening on is the interpreted
# plot data from the BNL website, which are at 300 K. 

# Temperature values to Doppler-Broaden to
Temps = [600, 900, 1200]

# Define the objects for containing data for H1, O16, U_235, U_238 
H1 = Nuclide_Data('H1', 1.008, [1, 1, 0])
O16 = Nuclide_Data('O16', 15.995, [1, 1, 0])
U235 = Nuclide_Data('U235', 235.044, [1, 1, 1])
U238 = Nuclide_Data('U238', 238.051, [1, 1, 1])

# Doppler-Broaden H1
for i in Temps:
    EM = [0, 0, 0]
    XS = [0, 0, 0]
    for j in range(len(H1.B)):
        if H1.B[j]:
            E1 = H1.EM[300][j]
            XS = H1.XS[300][j]
        H1.XS[i] = Doppler(EM, E1, XS, 300, i, H1.M)
        H1.EM[i] = EM

print('Finished Broadening H1')
