""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
"""

# ============================================================================
# Import statements
# ============================================================================
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})

# ============================================================================
# Classes
# ============================================================================
class isotope:
    def __init__(self, M):
        self.M = M # Atomic Mass
        
        self.file_contents = ['ES', 'NG']
        self.base_temp = 0 # K
        
        self.ESXS = {}
        self.ESEM = {}
        self.NGXS = {}
        self.NGEM = {}
        return
    
    def load_plot_data(self, filename):
        """ Loads BNL interpreted plotting data
        """
        ES, XS = np.loadtxt(filename, delimiter=',', unpack=True, skiprows=1)
        return ES, XS
    
    def load_data(self, filename, file_content, temperature):
        EM, XS = np.loadtxt(filename, delimiter=',', unpack=True, skiprows=1)
        if file_content == 'ES':
            self.ESXS[str(temperature)] = XS
            self.ESEM[str(temperature)] = EM
        elif file_content == 'NG':
            self.NGXS[str(temperature)] = XS
            self.NGEM[str(temperature)] = EM
            
# ============================================================================
# Load interpreted plotted data files from BNL website
# ============================================================================

# Plot H1 elastic scattering cross section library and Doppler Broadened
H1 = isotope(M=1.008)
H1.load_data('Data/H1_ES.txt', 'ES', 300)
H1.load_data('H1_ES_600', 'ES', 600)
H1.load_data('H1_ES_900', 'ES', 900)
H1.load_data('H1_ES_1200', 'ES', 1200)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
for key in H1.ESEM.keys():
    x = H1.ESEM[key]
    y = H1.ESXS[key]
    ax.loglog(x, y)
plt.xlabel('Energy [eV]'), plt.ylabel('Cross Section [Barns]')
legend = [i+' [K]' for i in list(H1.ESEM.keys())]
plt.legend(legend, loc=3), plt.title('H1 Elastic Scatter'), plt.show()

# Plot O16 elastic scattering cross section library and Doppler Broadened
O16 = isotope(M=15.995)
O16.load_data('Data/O16_ES.txt', 'ES', 300)
O16.load_data('O16_ES_600', 'ES', 600)
O16.load_data('O16_ES_900', 'ES', 900)
O16.load_data('O16_ES_1200', 'ES', 1200)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
for key in O16.ESEM.keys():
    x = O16.ESEM[key]
    y = O16.ESXS[key]
    ax.loglog(x, y)
legend = [i+' [K]' for i in list(O16.ESEM.keys())]
plt.xlabel('Energy [eV]'), plt.ylabel('Cross Section [Barns]')
plt.legend(legend, loc=3), plt.title('O16 Elastic Scatter'), plt.show()

# Plot U238 elastic scattering cross section library and Doppler Broadened
U238 = isotope(M=238.051)
U238.load_data('Data/U238_ES.txt', 'ES', 300)
U238.load_data('U238_ES_600', 'ES', 600)
U238.load_data('U238_ES_900', 'ES', 900)
U238.load_data('U238_ES_1200', 'ES', 1200)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
for key in U238.ESEM.keys():
    x = U238.ESEM[key]
    y = U238.ESXS[key]
    ax.loglog(x, y)
legend = [i+' [K]' for i in list(U238.ESEM.keys())]
plt.xlabel('Energy [eV]'), plt.ylabel('Cross Section [Barns]')
plt.legend(legend, loc=3), plt.title('U238  Elastic Scatter'), plt.show()

# Plot U238 radiative capture cross section library and Doppler Broadened
U238.load_data('Data/U238_NG.txt', 'NG', 300)
U238.load_data('U238_NG_600', 'NG', 600)
U238.load_data('U238_NG_900', 'NG', 900)
U238.load_data('U238_NG_1200', 'NG', 1200)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
for key in U238.NGEM.keys():
    x = U238.NGEM[key]
    y = U238.NGXS[key]
    ax.loglog(x, y)
legend = [i+' [K]' for i in list(U238.NGEM.keys())]
plt.xlabel('Energy [eV]'), plt.ylabel('Cross Section [Barns]')
plt.legend(legend, loc=3), plt.title('U238 Radiative Capture'), plt.show()

# Plot U235 elastic scattering cross section library and Doppler Broadened
U235 = isotope(M=235.044)
U235.load_data('Data/U235_ES.txt', 'ES', 300)
U235.load_data('U235_ES_600', 'ES', 600)
U235.load_data('U235_ES_900', 'ES', 900)
U235.load_data('U235_ES_1200', 'ES', 1200)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
for key in U235.ESEM.keys():
    x = U235.ESEM[key]
    y = U235.ESXS[key]
    ax.loglog(x, y)
legend = [i+' [K]' for i in list(U235.ESEM.keys())]
plt.xlabel('Energy [eV]'), plt.ylabel('Cross Section [Barns]')
plt.legend(legend, loc=3), plt.title('U235 Elastic Scatter'), plt.show()

# Plot U235 radiative Capture cross section library and Doppler Broadened
U235 = isotope(M=235.044)
U235.load_data('Data/U235_NG.txt', 'NG', 300)
U235.load_data('U235_NG_600', 'NG', 600)
U235.load_data('U235_NG_900', 'NG', 900)
U235.load_data('U235_NG_1200', 'NG', 1200)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
for key in U235.NGEM.keys():
    x = U235.NGEM[key]
    y = U235.NGXS[key]
    ax.loglog(x, y)
legend = [i+' [K]' for i in list(U235.NGEM.keys())]
plt.xlabel('Energy [eV]'), plt.ylabel('Cross Section [Barns]')
plt.legend(legend, loc=3), plt.title('U235 Radiative Capture'), plt.show()

# Compare BNL cross section data to BW DB cross section at 300
Emesh, y1 = np.loadtxt('Data/U238_NG_300_BW.txt', unpack=True, delimiter=',')
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
ax.loglog(Emesh, y1, U238.NGEM['300'], U238.NGXS['300']), plt.xlim(5, 12)
plt.legend(['$Breit-Wigner\ \sigma_\gamma$', '$BNL\ Interpreted\ Plot\ \sigma_\gamma$'], loc=3)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    