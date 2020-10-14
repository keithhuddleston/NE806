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
H1 = isotope(M=1.008)
H1.load_data('Data/H1_ES.txt', 'ES', 300)
H1.load_data('H1_ES_600', 'ES', 600)
H1.load_data('H1_ES_900', 'ES', 900)
H1.load_data('H1_ES_1200', 'ES', 1200)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
for key in H1.ESEM.keys():
    x = H1.ESEM[key]
    y = H1.ESXS[key]
    ax.loglog(x, y)
plt.xlabel('Energy [eV]'), plt.ylabel('Cross Section [Barns]')
legend = [i+' [K]' for i in list(H1.ESEM.keys())]
plt.legend(legend), plt.show()

O16 = isotope(M=15.995)
O16.load_data('Data/O16_ES.txt', 'ES', 300)
O16.load_data('O16_ES_600', 'ES', 600)
O16.load_data('O16_ES_900', 'ES', 900)
O16.load_data('O16_ES_1200', 'ES', 1200)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
for key in O16.ESEM.keys():
    x = O16.ESEM[key]
    y = O16.ESXS[key]
    ax.loglog(x, y)
legend = [i+' [K]' for i in list(O16.ESEM.keys())]
plt.legend(legend), plt.show()

U238 = isotope(M=238.051)
U238.load_data('Data/U238_ES.txt', 'ES', 300)
U238.load_data('U238_ES_600', 'ES', 600)
U238.load_data('U238_ES_900', 'ES', 900)
U238.load_data('U238_ES_1200', 'ES', 1200)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
for key in U238.ESEM.keys():
    x = U238.ESEM[key]
    y = U238.ESXS[key]
    ax.loglog(x, y)
legend = [i+' [K]' for i in list(U238.ESEM.keys())]
ax.set_xlim(1, 12)
plt.legend(legend), plt.show()
U238.load_data('Data/U238_NG.txt', 'NG', 300)
    
E, sigma_g = np.loadtxt('Data/U238_ES_0_BW.txt', unpack=True, delimiter=',')
Emesh, y1 = np.loadtxt('Data/U238_ES_300_BW.txt', unpack=True, delimiter=',')
Emesh, y2 = np.loadtxt('Data/U238_ES_600_BW.txt', unpack=True, delimiter=',')
Emesh, y3 = np.loadtxt('Data/U238_ES_900_BW.txt', unpack=True, delimiter=',')
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
ax.loglog(E, sigma_g, Emesh, y1, Emesh, y2, Emesh, y3), ax.set_xlim(1, 1000)
plt.legend(['$Breit-Wigner\ \sigma_\gamma, T_2=0K$', '$T_2 = 300K$', '$T_2 = 600K$', '$T_2 = 900K$'])
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    