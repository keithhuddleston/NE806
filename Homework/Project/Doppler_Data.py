""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
     
    Note, expect this file to run for awhile!
"""

# ============================================================================
# Import statements
# ============================================================================
import numpy as np
from Plot_BNL_Data import Nuclide_Data
from scipy.special import erf

# ============================================================================
# Create Text Files Containing the Doppler Broadened Cross Sections
# ============================================================================
def Doppler(E2, E1, S1, T1, T2, M, m=1.009):
    """ Doppler Broadens 

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
    M : Float
        Atomic mass of target.
    m : Float
        Atomic mass of projectile, 1.009 for Neutron.
    Returns
    -------
    S2 : array_like
        Reavaluated cross section data for energies E2, and temperature T2

    """
    Bk     = 6.617*10**-5        # Boltzman Constant, [eV K^-1]
    alpha  = (M/m)/(Bk*(T2-T1))  # Alpha term found in Doppler broadening Eqs.
    S2     = np.zeros(len(E2))   # Initialize new cross section data array
    S2_pos = np.zeros(len(E2))
    S2_neg = np.zeros(len(E2))
    
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
    
    def Af(E1, E2, S1, S2):
        den = (E2 - E1)
        num = (E2*S1) - (E1*S2)
        return num/den
    
    def Cf(E1, E2, S1, S2, alpha):
        den = (E2 - E1)*alpha
        num = (S2 - S1)
        return num/den
    
    # Evaluate Doppler-broadened cross section at specified energy E2[i]
    for i in range(len(E2)):
        S2i = 0
        y = [-1*np.sqrt(alpha*E2[i]), np.sqrt(alpha*E2[i])]
        for j in range(len(y)):
            Ek1 = E1[:-1]
            Ek2 = E1[1:]
            Sk1 = S1[:-1]
            Sk2 = S1[1:]            
            xk1 = np.sqrt(alpha*Ek1)
            xk2 = np.sqrt(alpha*Ek2)

            Ak = Af(Ek1, Ek2, Sk1, Sk2)
            Ck = Cf(Ek1, Ek2, Sk1, Sk2, alpha)

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

def Make_Doppler_Data(Nuclide, Temp, Energy=None):
    print('='*28 + '\nDoppler Broadening ' + Nuclide.N + ' Data' + '\n' + \
          '='*28 + '\n')
    data_type = ['ES', 'Total', 'Fission']
    for i in Temp:
        print('Broadening Data for Temperature ' + str(i) + \
              '[\N{DEGREE SIGN}K]...\n')
        NEM = [0, 0, 0]
        NXS = [0, 0, 0]
        for j in range(3):
            if Nuclide.B[j]:
                E1 = Nuclide.EM[300][j]
                XS = Nuclide.XS[300][j]
                NXS[j] = Doppler(E1, E1, XS, 300, i, Nuclide.M)
                NEM[j] = E1
                
        Nuclide.EM[i] = NEM
        Nuclide.XS[i] = NXS
        
        print('Saving Broadened Data to .txt files\n')
        for j in range(3):
            data = np.vstack(np.transpose([NEM[j], NXS[j]]))
            if Nuclide.B[j]:
                file_name = 'Data/Doppler/' + Nuclide.N + '_' + \
                    data_type[j] + '_' + str(i) + '.txt'
                np.savetxt(file_name, data, delimiter=',')
    print('Finished Broadening ' + Nuclide.N + ' Data\n')
    
# ============================================================================
# Load and Doppler-Broaden Data
# ============================================================================
# Note, the data we are performing Doppler Broadening on are the interpreted
# plot data from the BNL website, which are at 300 K. 
if __name__ == "__main__":
    # Temperature values to Doppler-Broaden to
    Temps = [600, 900, 1200]
    
    # Define the objects for containing data for H1, O16, U_235, U_238 
    H1 = Nuclide_Data('H1', 1.008, [1, 1, 0])
    O16 = Nuclide_Data('O16', 15.995, [1, 1, 0])
    U235 = Nuclide_Data('U235', 235.044, [1, 1, 1])
    U238 = Nuclide_Data('U238', 238.051, [1, 1, 1])
    
    # Doppler-Broaden Data
    Make_Doppler_Data(H1, Temps)
    Make_Doppler_Data(O16, Temps)
    Make_Doppler_Data(U235, Temps)
    Make_Doppler_Data(U238, Temps)

