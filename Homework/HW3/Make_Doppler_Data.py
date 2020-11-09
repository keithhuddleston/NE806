""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
     
    Note, expect this file to run for awhile!
"""

# ============================================================================
# Import statements
# ============================================================================
import numpy as np
from NE806_Functions import Doppler, Nuclide_Data # File written for class


# ============================================================================
# Create Text Files Containing the Doppler Broadened Cross Sections
# ============================================================================

# Note, the data we are performing Doppler Broadening on is the interpreted
# plot data from the BNL website, which are at 300 K. 

# Temperature values to Doppler-Broaden to
Temps = [600, 900, 1200]

# Energy mesh based on minimum and maximum values of the Casmo group structure
Emesh = np.logspace(-5, 7, 75000)

# Define the objects for containing data for H1, O16, U_235, U_238 
H1 = Nuclide_Data(1.008, False)
O16 = Nuclide_Data(15.995, False)
U235 = Nuclide_Data(235.044, True)
U238 = Nuclide_Data(238.051, True)

# Load in BNL plotted interpreted data for H1, O16, U_235, U_238 
H1.load_data('Data/BNL/H1_ES.txt', 'ES', 300)      # Elastic Scattering
O16.load_data('Data/BNL/O16_ES.txt', 'ES', 300)    # Elastic Scattering
U235.load_data('Data/BNL/U235_ES.txt', 'ES', 300)  # Elastic Scattering
U235.load_data('Data/BNL/U235_NG.txt', 'NG', 300)  # Radiative Capture
U238.load_data('Data/BNL/U238_ES.txt', 'ES', 300)  # Elastic Scattering
U238.load_data('Data/BNL/U238_NG.txt', 'NG', 300)  # Radiative Caputre

# Doppler-Broaden H1
for i in Temps:
    E1 = H1.ESEM['300']
    XS = H1.ESXS['300']
    H1.ESXS[str(i)] = Doppler(Emesh, E1, XS, 300, i, H1.M)
    H1.ESEM[str(i)] = Emesh
for key in list(H1.ESEM.keys())[1:]:
    np.savetxt('Data/Doppler/H1_ES_'+key, 
                np.vstack(np.transpose([H1.ESEM[key], 
                                        H1.ESXS[key]])), 
                                        delimiter=',')
print('Finished Broadening H1')

# Doppler-Broaden O16
for i in Temps:
    E1 = O16.ESEM['300']
    XS = O16.ESXS['300']
    O16.ESXS[str(i)] = Doppler(Emesh, E1, XS, 300, i, O16.M)
    O16.ESEM[str(i)] = Emesh
for key in list(O16.ESEM.keys())[1:]:
    np.savetxt('Data/Doppler/O16_ES_'+key, 
                np.vstack(np.transpose([O16.ESEM[key], 
                                        O16.ESXS[key]])), 
                                        delimiter=',')
print('Finished Broadening O16')
  
# Doppler-Broaden U238
for i in Temps:
    E1 = U238.ESEM['300']
    XS = U238.ESXS['300']
    U238.ESXS[str(i)] = Doppler(Emesh, E1, XS, 300, i, U238.M)
    U238.ESEM[str(i)] = Emesh
    for key in list(U238.ESEM.keys())[1:]:
        np.savetxt('Data/Doppler/U238_ES_'+key, 
                    np.vstack(np.transpose([U238.ESEM[key], 
                                            U238.ESXS[key]])), 
                                            delimiter=',')
    E1 = U238.NGEM['300']
    XS = U238.NGXS['300']
    U238.NGXS[str(i)] = Doppler(Emesh, E1, XS, 300, i, U238.M)
    U238.NGEM[str(i)] = Emesh
    for key in list(U238.NGEM.keys())[1:]:
        np.savetxt('Data/Doppler/U238_NG_'+key, 
                    np.vstack(np.transpose([U238.NGEM[key], 
                                            U238.NGXS[key]])), 
                                            delimiter=',')
print('Finished Broadening U238')

# Doppler-Broaden U235
for i in Temps:
    E1 = U235.ESEM['300']
    XS = U235.ESXS['300']
    U235.ESXS[str(i)] = Doppler(Emesh, E1, XS, 300, i, U235.M)
    U235.ESEM[str(i)] = Emesh
    for key in list(U235.ESEM.keys())[1:]:
        np.savetxt('Data/Doppler/U235_ES_'+key, 
                    np.vstack(np.transpose([U235.ESEM[key], 
                                            U235.ESXS[key]])), 
                                            delimiter=',')        
    E1 = U235.NGEM['300']
    XS = U235.NGXS['300']
    U235.NGXS[str(i)] = Doppler(Emesh, E1, XS, 300, i, U235.M)
    U235.NGEM[str(i)] = Emesh
    for key in list(U235.NGEM.keys())[1:]:
        np.savetxt('Data/Doppler/U235_NG_'+key, 
                    np.vstack(np.transpose([U235.NGEM[key], 
                                            U235.NGXS[key]])), 
                                            delimiter=',') 
print('Finished Broadening U235')

# ============================================================================
# Get Original BNL Data On Same Energy Mesh
# ============================================================================
H1.ESXS['300'] = np.interp(Emesh, H1.ESEM['300'], H1.ESXS['300'])
H1.ESEM['300'] = Emesh
np.savetxt('Data/Doppler/H1_ES_300.txt', 
           np.vstack(np.transpose([H1.ESEM['300'], 
                                   H1.ESXS['300']])), 
                                   delimiter=',') 

O16.ESXS['300'] = np.interp(Emesh, O16.ESEM['300'], O16.ESXS['300'])
O16.ESEM['300'] = Emesh
np.savetxt('Data/Doppler/O16_ES_300.txt', 
           np.vstack(np.transpose([O16.ESEM['300'], 
                                   O16.ESXS['300']])), 
                                   delimiter=',') 

U238.ESXS['300'] = np.interp(Emesh, U238.ESEM['300'], U238.ESXS['300'])
U238.ESEM['300'] = Emesh
U238.NGXS['300'] = np.interp(Emesh, U238.NGEM['300'], U238.NGXS['300'])
U238.NGEM['300'] = Emesh
np.savetxt('Data/Doppler/U238_ES_300.txt', 
           np.vstack(np.transpose([U238.ESEM['300'], 
                                   U238.ESXS['300']])), 
                                   delimiter=',') 
np.savetxt('Data/Doppler/U238_NG_300.txt', 
           np.vstack(np.transpose([U238.NGEM['300'], 
                                   U238.NGXS['300']])), 
                                   delimiter=',')

U235.ESXS['300'] = np.interp(Emesh, U235.ESEM['300'], U235.ESXS['300'])
U235.ESEM['300'] = Emesh
U235.NGXS['300'] = np.interp(Emesh, U235.NGEM['300'], U235.NGXS['300'])
U235.NGEM['300'] = Emesh
np.savetxt('Data/Doppler/U235_ES_300.txt', 
           np.vstack(np.transpose([U235.ESEM['300'], 
                                   U235.ESXS['300']])), 
                                   delimiter=',') 
np.savetxt('Data/Doppler/U235_NG_300.txt', 
           np.vstack(np.transpose([U235.NGEM['300'], 
                                   U235.NGXS['300']])), 
                                   delimiter=',') 