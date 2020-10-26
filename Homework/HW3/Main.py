""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
"""

# ============================================================================
# Import statements
# ============================================================================
import numpy as np
import matplotlib.pyplot as plt
from NE806_Functions import doppler # File written for class


# ============================================================================
# Classes
# ============================================================================
def unique(sequence):
    seen = set()
    return [x for x in sequence if not (x in seen or seen.add(x))]

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
    
    def load_data(self, filename, file_content, temperature):
        EM, XS = np.loadtxt(filename, delimiter=',', unpack=True, skiprows=1)
        if len(EM) != len(set(EM)):
            ind = np.zeros(len(EM)-len(set(EM)), dtype=int)
            if len(ind) > 0:
                j = 0
                for i in range(len(EM[:-1])):
                    if EM[i] == EM[i+1]:
                        ind[j] = i
                        j += 1
            EM = np.delete(EM, ind)
            XS = np.delete(XS, ind)
        if file_content == 'ES':
            self.ESXS[str(temperature)] = XS
            self.ESEM[str(temperature)] = EM
        elif file_content == 'NG':
            self.NGXS[str(temperature)] = XS
            self.NGEM[str(temperature)] = EM


# ============================================================================
# Functions
# ============================================================================

# Narrow Resonance Flux Approximation
def phinr(sigma_t, sigma_o, E):
    return 1 / (E*(sigma_t + sigma_o))

# Wide Resonance Flux Approximation
def phiwr(sigma_a, sigma_o, E):
    return 1 / (E*(sigma_a + sigma_o))

# Microscopic Flux Weighted Average Cross Section
def sigma_g(Phi, XS, Emesh, Group_Structure):
    return(sum(Phi*XS)/sum(Phi))

# def sigma_g(Phi, XS, Emesh, Group_Structure):
#     return(np.trapz(Phi*XS)/np.trapz(Phi))

def make_group(Emesh, phi, xs, group_structure):
    ind_G = np.zeros(len(group_structure), dtype=int)
    k = 0
    for i in range(len(Emesh)):
        if Emesh[i] > c16[k]:
            ind_G[k] = int(i)
            k += 1
    if ind_G[-1] == 0:
        ind_G[-1] = len(Emesh)

    group_e = [] # Group-wise energy values
    group_s = [] # Group-wise cross section values
    group_p = [] # Group-wise flux values
    
    c16_s = np.interp(c16, Emesh, xs)  # Interpolate boundary values
    c16_p = np.interp(c16, Emesh, phi) # Interpolate boundary values
    
    for i in range(len(ind_G[:-1])):
        L1 = [c16[i]] + list(Emesh[ind_G[i]:ind_G[i+1]]) + [c16[i+1]]
        L2 = [c16_s[i]] + list(xs[ind_G[i]:ind_G[i+1]]) + [c16_s[i+1]]
        L3 = [c16_p[i]] + list(phi[ind_G[i]:ind_G[i+1]]) + [c16_p[i+1]]
        group_e.append(np.array(L1))
        group_s.append(np.array(L2))
        group_p.append(np.array(L3))
    if group_e[-1][-1] == group_e[-1][-2]:
        group_e[-1] = group_e[-1][:-1]
        group_s[-1] = group_s[-1][:-1]
        group_p[-1] = group_p[-1][:-1]
    return group_e, group_s, group_p


# ============================================================================
# Create Doppler Broadened Cross Sections
# ============================================================================
# Note this takes forever to run, was used to make doppler broadened xs
Flag = False
if Flag:
    Temps = [600, 900, 1200] # Temperature values to Doppler-Broaden at
    
    # Emesh based on Casmo group structure
    Emesh = np.logspace(-5, 7, 75000)
    
    # Load in BNL plotted interpreted data at for H1
    # Elastic Scattering
    H1 = isotope(M=1.008)
    H1.load_data('Data/H1_ES.txt', 'ES', 300)
    
    # Load in BNL plotted interpreted data at for O16
    # Elastic Scattering
    O16 = isotope(M=15.995)
    O16.load_data('Data/O16_ES.txt', 'ES', 300)
    
    # Load in BNL plotted interpreted data at for U235
    # Elastic Scattering
    # Radiative Capture
    U235 = isotope(M=235.044)
    U235.load_data('Data/U235_ES.txt', 'ES', 300)
    U235.load_data('Data/U235_NG.txt', 'NG', 300)
    
    # Load in BNL plotted interpreted data at for U238
    # Elastic Scattering
    # Radiative Capture
    U238 = isotope(M=238.051)
    U238.load_data('Data/U238_ES.txt', 'ES', 300)
    U238.load_data('Data/U238_NG.txt', 'NG', 300)
    
    # Doppler broaden all loaded data
    # Doppler Broaden H1
    for i in Temps:
        E1 = H1.ESEM['300']
        XS = H1.ESXS['300']
        H1.ESXS[str(i)] = doppler(Emesh, E1, XS, 300, i, H1.M)
        H1.ESEM[str(i)] = Emesh
    for key in list(H1.ESEM.keys())[1:]:
        np.savetxt('H1_ES_'+key, np.vstack(np.transpose([H1.ESEM[key], 
                                                         H1.ESXS[key]])), 
                                                         delimiter=',')
    
    # Doppler Broaden O16
    for i in Temps:
        E1 = O16.ESEM['300']
        XS = O16.ESXS['300']
        O16.ESXS[str(i)] = doppler(Emesh, E1, XS, 300, i, O16.M)
        O16.ESEM[str(i)] = Emesh
    for key in list(O16.ESEM.keys())[1:]:
        np.savetxt('O16_ES_'+key, np.vstack(np.transpose([O16.ESEM[key], 
                                                          O16.ESXS[key]])), 
                                                          delimiter=',')
       
    # Doppler Broaden U238
    for i in Temps:
        E1 = U238.ESEM['300']
        XS = U238.ESXS['300']
        U238.ESXS[str(i)] = doppler(Emesh, E1, XS, 300, i, U238.M)
        U238.ESEM[str(i)] = Emesh
        for key in list(U238.ESEM.keys())[1:]:
            np.savetxt('U238_ES_'+key, np.vstack(np.transpose([U238.ESEM[key], 
                                                               U238.ESXS[key]])), 
                                                               delimiter=',')
        E1 = U238.NGEM['300']
        XS = U238.NGXS['300']
        U238.NGXS[str(i)] = doppler(Emesh, E1, XS, 300, i, U238.M)
        U238.NGEM[str(i)] = Emesh
        for key in list(U238.NGEM.keys())[1:]:
            np.savetxt('U238_NG_'+key, np.vstack(np.transpose([U238.NGEM[key], 
                                                               U238.NGXS[key]])), 
                                                               delimiter=',')

    # Doppler Broaden U235
    for i in Temps:
        E1 = U235.ESEM['300']
        XS = U235.ESXS['300']
        U235.ESXS[str(i)] = doppler(Emesh, E1, XS, 300, i, U235.M)
        U235.ESEM[str(i)] = Emesh
        for key in list(U235.ESEM.keys())[1:]:
            np.savetxt('U235_ES_'+key, np.vstack(np.transpose([U235.ESEM[key], 
                                                               U235.ESXS[key]])), 
                                                               delimiter=',')        
        E1 = U235.NGEM['300']
        XS = U235.NGXS['300']
        U235.NGXS[str(i)] = doppler(Emesh, E1, XS, 300, i, U235.M)
        U235.NGEM[str(i)] = Emesh
        for key in list(U235.NGEM.keys())[1:]:
            np.savetxt('U235_NG_'+key, np.vstack(np.transpose([U235.NGEM[key], 
                                                               U235.NGXS[key]])), 
                                                               delimiter=',')       


# ============================================================================
# Testing, Litmus test
# ============================================================================
# Casmo 2 group structure
c2 = np.array([.001e1, 1.00e-7, 1.00e-11])*1e6 # eV

# Casmo 16 group structure
c16 = np.array([1.00e1 , 8.21e-1, 5.53e-3, 4.00e-6, 1.30e-6, 1.15e-6, 1.097e-6,
                1.02e-6, 9.71e-7, 8.50e-7, 6.25e-7, 3.50e-7, 2.80e-7, 1.40e-7 ,
                5.80e-8, 3.00e-8, 1.00e-11])*1e6 # eV

c2 = c2[::-1]
c16 = c16[::-1]

# Microscopic Dilution/Background Cross Section
sd = [1e1, 1e2, 1e3, 1e4, 1e5]

if __name__ == "__main__":
    e = np.logspace(-5, 7, 75000)
    p = np.ones(len(e))
    s = np.ones(len(e))
    e, s, p = make_group(e, p, s, c16)
    sg = np.zeros(len(e))
    for i in range(len(sg)):
        sg[i] = sigma_g(p[i], s[i], e[i], c16)
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
    ax.set_xscale( "log" )
    ax.set_yscale( "log" )
    for i in range(len(sg)):
        ax.fill_between(e[i], sg[i])


# ============================================================================
# Testing, Answer Homework Questions
# ============================================================================
Flag = True
if Flag:
    # Note, this block of text won't run unless Doppler Broadening was done
    
    # Load H1 Elastic Scattering Data
    H1 = isotope(M=1.008)
    H1.load_data('Data/H1_ES.txt', 'ES', 300)
    H1.load_data('H1_ES_600', 'ES', 600)
    H1.load_data('H1_ES_900', 'ES', 900)
    H1.load_data('H1_ES_1200', 'ES', 1200)
    
    # Load O16 Elastic Scattering Data
    O16 = isotope(M=15.995)
    O16.load_data('Data/O16_ES.txt', 'ES', 300)
    O16.load_data('O16_ES_600', 'ES', 600)
    O16.load_data('O16_ES_900', 'ES', 900)
    O16.load_data('O16_ES_1200', 'ES', 1200)
    
    # Load U238 Elastic Scattering Data
    U238 = isotope(M=238.051)
    U238.load_data('Data/U238_ES.txt', 'ES', 300)
    U238.load_data('U238_ES_600', 'ES', 600)
    U238.load_data('U238_ES_900', 'ES', 900)
    U238.load_data('U238_ES_1200', 'ES', 1200)
    
    # Load U238 Radiative Capture Data
    U238.load_data('Data/U238_NG.txt', 'NG', 300)
    U238.load_data('U238_NG_600', 'NG', 600)
    U238.load_data('U238_NG_900', 'NG', 900)
    U238.load_data('U238_NG_1200', 'NG', 1200)
    
    # Load U235 Elastic Scattering Data
    U235 = isotope(M=235.044)
    U235.load_data('Data/U235_ES.txt', 'ES', 300)
    U235.load_data('U235_ES_600', 'ES', 600)
    U235.load_data('U235_ES_900', 'ES', 900)
    U235.load_data('U235_ES_1200', 'ES', 1200)
    
    # Load U235 Radiative Capture
    U235.load_data('Data/U235_NG.txt', 'NG', 300)
    U235.load_data('U235_NG_600', 'NG', 600)
    U235.load_data('U235_NG_900', 'NG', 900)
    U235.load_data('U235_NG_1200', 'NG', 1200)
    
    Emesh   = U235.ESEM['600']
    sigma_a = U235.NGXS['600']
    xs      = U235.NGXS['600']
    sigma_e = U235.ESXS['600']
    phi     = phinr(sigma_a+sigma_e, sd[0], Emesh)
    e, s, p = make_group(Emesh, phi, xs, c16)
    
    sg = np.zeros(len(e))
    
    for i in range(len(sg)):
        sg[i] = sigma_g(p[i], s[i], e[i], c16)
        
    # Plot Results
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
    for i in range(len(e)):
        ax.loglog(e[i], p[i])
        ax.axvline(c16[i], ls='--', c='k')
    plt.plot()
    
    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
    ax1.set_xscale( "log" )
    ax1.set_yscale( "log" )
    # ax2 = ax1.twinx()
    for i in range(len(sg)):
        ax1.fill_between(e[i], sg[i])
        ax1.loglog(e[i], s[i], color='black', linestyle='--')
        # ax1.loglog(e[i], p[i], color='black')
    # ax1.loglog(e[-3], s[-1], color='black', linestyle='--', maker='o')

