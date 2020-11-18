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

def Doppler_Nuclide(Nuclide, Temp, Energy)
    for i in Temp:
        NEM = [0, 0, 0]
        NXS = [0, 0, 0]
        for j in range(len(Nuclide.B)):
            if H1.B[j]:
                E1 = Nuclide.EM[300][j]
                XS = Nuclide.XS[300][j]
                NXS[j] = Doppler(E1, E1, XS, 300, i, Nuclide.M)
                NEM[j] = E1
                
        Nuclide.EM[i] = NEM
        Nuclide.XS[i] = NXS
        for j in NEM:
            np.savetxt('Data/Doppler/' + H1.N +'_ES_'+key, 
                np.vstack(np.transpose([H1.ESEM[key], 
                                        H1.ESXS[key]])), 
                                        delimiter=',')

print('Finished Broadening H1')
