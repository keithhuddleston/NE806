""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
"""

# ============================================================================
# Import Statements
# ============================================================================
import numpy as np
import matplotlib.pyplot as plt
from NE806_Functions import Nuclide_Data
plt.rcParams.update({'font.size': 18})

            
# ============================================================================
# Load Interpolated Interpreted Plotted Data Files and Doppler-Broadened Data
# ============================================================================
H1 = Nuclide_Data(1.008, False)
O16 = Nuclide_Data(15.995, False)
U238 = Nuclide_Data(238.051, True)
U235 = Nuclide_Data(235.044, True)

H1.load_data('Data/Doppler/H1_ES_300.txt', 'ES', 300)
H1.load_data('Data/Doppler/H1_ES_600',     'ES', 600)
H1.load_data('Data/Doppler/H1_ES_900',     'ES', 900)
H1.load_data('Data/Doppler/H1_ES_1200',    'ES', 1200)

O16.load_data('Data/Doppler/O16_ES_300.txt', 'ES', 300)
O16.load_data('Data/Doppler/O16_ES_600',     'ES', 600)
O16.load_data('Data/Doppler/O16_ES_900',     'ES', 900)
O16.load_data('Data/Doppler/O16_ES_1200',    'ES', 1200)

U238.load_data('Data/Doppler/U238_ES_300.txt', 'ES', 300)
U238.load_data('Data/Doppler/U238_ES_600',     'ES', 600)
U238.load_data('Data/Doppler/U238_ES_900',     'ES', 900)
U238.load_data('Data/Doppler/U238_ES_1200',    'ES', 1200)
U238.load_data('Data/Doppler/U238_NG_300.txt', 'NG', 300)
U238.load_data('Data/Doppler/U238_NG_600',     'NG', 600)
U238.load_data('Data/Doppler/U238_NG_900',     'NG', 900)
U238.load_data('Data/Doppler/U238_NG_1200',    'NG', 1200)

U235.load_data('Data/Doppler/U235_ES_300.txt', 'ES', 300)
U235.load_data('Data/Doppler/U235_ES_600',     'ES', 600)
U235.load_data('Data/Doppler/U235_ES_900',     'ES', 900)
U235.load_data('Data/Doppler/U235_ES_1200',    'ES', 1200)
U235.load_data('Data/Doppler/U235_NG_300.txt', 'NG', 300)
U235.load_data('Data/Doppler/U235_NG_600',     'NG', 600)
U235.load_data('Data/Doppler/U235_NG_900',     'NG', 900)
U235.load_data('Data/Doppler/U235_NG_1200',    'NG', 1200)

# ============================================================================
# Plotting Time
# ============================================================================
# Plot H1 Elastic Scattering Data
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
for key in H1.ESEM.keys():
    x = H1.ESEM[key]
    y = H1.ESXS[key]
    ax.loglog(x, y)
plt.xlabel('Energy [eV]')
plt.ylabel('Cross Section [Barns]')
plt.title('H1 Elastic Scatter')
legend = [i+' [K]' for i in list(H1.ESEM.keys())]
plt.legend(legend, loc=3)
plt.show()

# Plot O16 Elastic Scattering Data
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
for key in O16.ESEM.keys():
    x = O16.ESEM[key]
    y = O16.ESXS[key]
    ax.loglog(x, y)
legend = [i+' [K]' for i in list(O16.ESEM.keys())]
plt.xlabel('Energy [eV]'), plt.ylabel('Cross Section [Barns]')
plt.legend(legend, loc=3), plt.title('O16 Elastic Scatter'), plt.show()

# Plot U238 Elastic Scattering Data
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
for key in U238.ESEM.keys():
    x = U238.ESEM[key]
    y = U238.ESXS[key]
    ax.loglog(x, y)
legend = [i+' [K]' for i in list(U238.ESEM.keys())]
plt.xlabel('Energy [eV]'), plt.ylabel('Cross Section [Barns]')
plt.legend(legend, loc=3), plt.title('U238  Elastic Scatter'), plt.show()

# Plot U238 Radiative Capture Data
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
for key in U238.NGEM.keys():
    x = U238.NGEM[key]
    y = U238.NGXS[key]
    ax.loglog(x, y)
legend = [i+' [K]' for i in list(U238.NGEM.keys())]
plt.xlabel('Energy [eV]'), plt.ylabel('Cross Section [Barns]')
plt.legend(legend, loc=3), plt.title('U238 Radiative Capture'), plt.show()

# Plot U235 Elastic Scattering Data
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
for key in U235.ESEM.keys():
    x = U235.ESEM[key]
    y = U235.ESXS[key]
    ax.loglog(x, y)
legend = [i+' [K]' for i in list(U235.ESEM.keys())]
plt.xlabel('Energy [eV]'), plt.ylabel('Cross Section [Barns]')
plt.legend(legend, loc=3), plt.title('U235 Elastic Scatter'), plt.show()

# Plot U235 Radiative Capture Data
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
for key in U235.NGEM.keys():
    x = U235.NGEM[key]
    y = U235.NGXS[key]
    ax.loglog(x, y)
legend = [i+' [K]' for i in list(U235.NGEM.keys())]
plt.xlabel('Energy [eV]'), plt.ylabel('Cross Section [Barns]')
plt.legend(legend, loc=3), plt.title('U235 Radiative Capture'), plt.show()

# Compare BNL Cross Section Data to BW DB cross section at 300
Emesh, y1 = np.loadtxt('Data/Doppler/U238_NG_300_BW.txt', 
                       unpack=True, delimiter=',')
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
ax.loglog(Emesh, y1, U238.NGEM['300'], U238.NGXS['300']), plt.xlim(5, 12)
plt.legend(['$Breit-Wigner\ \sigma_\gamma$', 
            '$BNL\ Interpreted\ Plot\ \sigma_\gamma$'], loc=3)
