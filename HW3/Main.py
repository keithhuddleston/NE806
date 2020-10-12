""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
"""

# ============================================================================
# Import statements
# ============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
import re

# ============================================================================
# Doppler Broadening Equations
# ============================================================================
F0 = lambda a: erf(a)
H0 = lambda a, b: F0(a) - F0(b)

F1 = lambda a: np.sqrt(1/np.pi) * (1-np.exp(-a**2))
H1 = lambda a, b: F1(a) - F1(b)

F2 = lambda a: (1/2)*erf(a) - (a/np.sqrt(np.pi))*np.exp(-a**2)
H2 = lambda a, b: F2(a) - F2(b)

F3 = lambda a: np.sqrt(1/np.pi) * (1-(1+a**2)*np.exp(-a**2))
H3 = lambda a, b: F3(a) - F3(b)

F4 = lambda a: (3/4)*erf(a) - np.sqrt(1/np.pi)*((3*a/2)+a**3)*np.exp(-a**2)
H4 = lambda a, b: F4(a) - F4(b)

def Akf(E1, E2, S1, S2):
    den = (E2 - E1)
    num = (E2*S1) - (E1*S2)
    return num/den

def Ckf(E1, E2, S1, S2, alpha):
    den = (E2 - E1)*alpha
    num = (S2 - S1)
    return num/den

def doppler(E2, E1, S1, T1, T2, M, m=1.009,  neg=False):
    """

    Parameters
    ----------
    E2 : array_like
        New energy mesh at which to evaluate cross sections at.
    E1 : array_like
        Energy mesh of the reference cross section data.
    S1 : array_like
        Cross section data of reference.
    T1 : Float
        Temperature of the reference cross section data.
    T2 : Float
        New temperature at which to evaluate cross sections at.
    m : Float
        Atomic mass of projectile, 1.009 for Neutron.
    M : Float
        Atomic mass of target.
    neg : Int, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    S2 : array_like
        Reavaluated cross section data for energies E2, and temperature T2

    """
    Bk = 6.617*10**-5           # Boltzman Constant, [ev K^-1]
    alpha = (M/m)/(Bk*(T2-T1))  # Alpha term found in Doppler broadening Eqs.
    S2 = np.zeros(len(E2))      # Initialize new cross section data array
    S2_pos = np.zeros(len(E2)) 
    S2_neg = np.zeros(len(E2)) 
    # Evaluate Doppler-broadened cross section at specified energy E2[i]
    for i in range(len(E2)):
        S2i = 0

        E2i = E2[i]
        y = np.sqrt(alpha*E2i)
        y = [-1*y, 1*y]
        for j in range(len(y)):
            Ek1 = E1[:-1]
            Ek2 = E1[1:]
            
            xk1 = np.sqrt(alpha*Ek1)
            xk2 = np.sqrt(alpha*Ek2)
    
            Sk1 = S1[:-1]
            Sk2 = S1[1:]
    
            Ak = Akf(Ek1, Ek2, Sk1, Sk2)
            Ck = Ckf(Ek1, Ek2, Sk1, Sk2, alpha)
    
            Zk1 = xk1 - y[j]
            Zk2 = xk2 - y[j]
    
            H0k = H0(Zk2, Zk1)
            H1k = H1(Zk2, Zk1)
            H2k = H2(Zk2, Zk1)
            H3k = H3(Zk2, Zk1)
            H4k = H4(Zk2, Zk1)
    
            S2i = H4k * (Ck) \
                 + H3k * (4*Ck*y[j]) \
                 + H2k * (Ak+6*Ck*y[j]**2) \
                 + H1k * (2*Ak*y[j]+4*Ck*y[j]**3) \
                 + H0k * (Ak*y[j]**2+Ck*y[j]**4)
            S2i = sum(S2i)
            if j == 0:
                S2_neg[i] = S2i/2/y[j]**2
            else:
                S2_pos[i] = S2i/2/y[j]**2
            S2 = S2_pos - S2_neg
    return S2

# ============================================================================
# Load in BNL interpreted data files from the Data subfolder
# ============================================================================
# Casmo 2 group structure
c2 = np.array([.001e1, 1.00e-7, 1.00e-11])*1e6 # eV

# Casmo 16 group structure
c16 = np.array([1.00e1 , 8.21e-1, 5.53e-3, 4.00e-6, 1.30e-6, 1.15e-6, 1.097e-6,
                1.02e-6, 9.71e-7, 8.50e-7, 6.25e-7, 3.50e-7, 2.80e-7, 1.40e-7 ,
                5.80e-8, 3.00e-8, 1.00e-11])*1e6 # eV

def load_BNL_RP(filename):
    """ File for loading BNL interpreted resonance parameter files
    """
    file = open(filename, 'r')
    lines = file.readlines()
    l = len(lines)
    data = np.zeros((l, 6)) # There should always be six columns in this file
    
    # Alter file values so that python can recognize them as floating point
    p1 = r'(\+*\-*\d\.\d*\+*\-*\d+)'
    for i in range(l):
        line = lines[i]
        found = re.findall(p1, line)
        for j in range(6):
            value = found[j]
            for k in range(len(value))[::-1]:
                if value[k] == '-' or value[k] == '+':
                    value = value[:k]+'e'+value[k:]
                    break
            data[i][j] = value
            
    # Seperate data into column vectors containing E0, J, GN, GG, GFA, GFB
    E0, J, GN, GG, GFA, GFB = [np.zeros(l) for i in range(6)]
    for i in range(l):
        E0[i], J[i], GN[i], GG[i], GFA[i], GFB[i] = data[i]
        ind = 0
    for i in range(len(E0)):
        if E0[i] < 0:
            ind = i+1
    return E0[ind:], J[ind:], GN[ind:], GG[ind:], GFA[ind:], GFB[ind:]

def breitWigner(E0, J, GN, GG, GFA, GFB, A, sigma_p, I, E):
    """ Reconstructs cross section data given resonance parameters
    """
    G = GN+GG
    sigma_g = np.zeros_like(E)
    sigma_n = np.zeros_like(E)
    sigma_n[:] = sigma_p
    sigma_e = np.zeros_like(E)
    #  Nuclear Radius in (cm), units cancel out with scattering interfence term
    R = 1.25*10**-13 * 238**(1/3)
    
    for i in range(len(E0)):
        # statistical spin factor, where I=nuclear spint, J=total spin
        g = (2*J[i]+1)/(2*(2*I+1))
    
        # total cross section at resonance energy (DH 2-36)
        sigma_0 = 2.608e6 * (A+1)**2/(A**2 * E0[i]) * (GN[i]/G[i]) * g
    
        # capture cross section (DH 2-35)
        y = (2/G[i])*(E - E0[i])
        sigma_g += sigma_0 * (GG[i]/G[i]) * np.sqrt(E0[i]/E) * (1/(1+y**2))
        sigma_e += sigma_0 * (GN[i]/G[i]) * np.sqrt(E0[i]/E) * (1/(1+y**2)) \
                 + sigma_0 * (2*R/(4.55*10**-10/E0[i]**0.5))*(y/(1+y**2))

    sigma_e = sigma_e + sigma_n  # Add potential scattering
    return sigma_e, sigma_g

def fwxs(E, X, Phi, GS):
    """
    Input:
          E, energy mesh of X and Phi
          X, Cross section data
          Phi, flux
          GS, specified group structure
    Output:
          GX,
    """
    # assert len(E) == len(X) == len(Phi), 'Length of input must be the same'
    # Step 1 add e values of GS to E
    # Step 2 interpolate values for new Es
    # Step 3 perform integral
    plt.loglog(E, X)
    plt.vlines(GS, min(X), max(X), linestyle='--'), plt.show()

# ============================================================================
# Classes
# ============================================================================
class isotope:
    def __init__(self, M):
        self.M = M # Atomic Mass
        
        self.file_formats = ['BNL_plot', 'BNL_resonance']
        self.file_contents = ['ES', 'NG']
        self.tempbase = 0 # K
        
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
    
    def load_data(self, filename, file_format, file_content):
        if file_format == 'BNL_plot':
            EM, XS = self.load_plot_data(filename)
            self.tempbase = 300 # K
        else:
            print('Not ready to load resonance parameter files :(')
            self.tempbase = 0 # K
        if file_content == 'ES':
            self.ESXS[str(self.tempbase)] = XS
            self.ESEM[str(self.tempbase)] = EM
        elif file_content == 'NG':
            self.NGXS[str(self.tempbase)] = XS
            self.NGEM[str(self.tempbase)] = EM
      
# ============================================================================
# Testing Oct 9. 2020
# ============================================================================
Flag = False
if Flag:
    # Test 1. Check that function form of BW is working right
    E, XS = np.loadtxt('Data/U238_NG.txt', delimiter=',', unpack=True, skiprows=1)
    
    I = 0 
    A = 238
    sigma_p = 11.2934
    E1 = np.logspace(-11, 10, 50000)
    E0, J, GN, GG, GFA, GFB = load_BNL_RP('Data/U238_RP_L0.txt')
    sigma_e, sigma_g = breitWigner(E0, J, GN, GG, GFA, GFB, A, sigma_p, I, E1)
    
    # Compare Dr. Robert's resonance file to full file
    E01, J1, GN1, GG1, GFA1, GFB1 = np.loadtxt('u238.txt', unpack=True)
    
    # Resonance recontruction using reduced file
    sigma_e0, sigma_g0 = breitWigner(E01, J1, GN1, GG1, GFA1, 
                                     GFB1, A, sigma_p, I, E1)
    
    # Do again for l=1 resonance data
    E0, J, GN, GG, GFA, GFB = load_BNL_RP('Data/U238_RP_L1.txt')
    sigma_e1, sigma_g1 = breitWigner(E0, J, GN, GG, GFA, GFB, A, sigma_p, I, E1)
    
    # Plot results
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(12,24))
    ax1.loglog(E, XS, E1, sigma_g0, E1, sigma_g, E1, sigma_g1)
    ax1.legend(['BNL interpreted', 'Reduced l=0', 'Full l=0', 'Full l=1'])
    ax1.set_xlim(5, 10)
    ax1.set_ylabel('Cross Section (Barns)')
    
    ax2.loglog(E, XS, E1, sigma_g0, E1, sigma_g, E1, sigma_g1)
    ax2.set_xlim(700, 1000)
    ax2.set_ylabel('Cross Section (Barns)')
    
    ax3.loglog(E, XS, E1, sigma_g0, E1, sigma_g, E1, sigma_g1)
    ax3.set_xlim(1000, 1300)
    ax3.set_ylabel('Cross Section (Barns)')
    
    ax4.loglog(E, XS, E1, sigma_g0, E1, sigma_g, E1, sigma_g1)
    ax4.vlines(c16[1:], ymin=min(sigma_e), ymax=max(sigma_e), linestyle='--')
    ax4.set_xlim(min(c16), max(c16))
    ax4.set_xlabel('Energy (eV)')
    ax4.set_ylabel('Cross Section (Barns)')
    plt.show()
    
# ============================================================================
# Testing Oct 12. 2020
# ============================================================================
Flag = True
if Flag:
    Temps = [0, 300, 600, 900, 1200] # Temperature values specified in HW
    
    H1i = isotope(M=1.008)
    H1i.load_data('Data/H1_ES.txt', 'BNL_plot', 'ES')
    
    O16 = isotope(M=15.995)
    O16.load_data('Data/O16_ES.txt', 'BNL_plot', 'ES')
    
    U235 = isotope(M=235.044)
    U235.load_data('Data/U235_ES.txt', 'BNL_plot', 'ES')
    U235.load_data('Data/U235_NG.txt', 'BNL_plot', 'NG')
    
    U238 = isotope(M=238.051)
    U238.load_data('Data/U238_ES.txt', 'BNL_plot', 'ES')
    U238.load_data('Data/U238_NG.txt', 'BNL_plot', 'NG')

m = 1.009
M = 1.008
Emesh = H1i.ESEM['300']
E = H1i.ESEM['300']
XS = H1i.ESXS['300']
H1i.ESXS['600'] = doppler(Emesh, E, XS, 300, 600, m, M)
H1i.ESEM['600'] = E
for key in H1i.ESEM.keys():
    print(key)
    plt.loglog(H1i.ESEM[key], H1i.ESXS[key])






























