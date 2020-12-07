""" Final Project for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
    
    Contains the user written classes or functions that are used in multiple 
    other files.
    
    Last Edit: Dec. 7, 2020
"""

# ============================================================================
# Import Statements
# ============================================================================
import numpy as np
import os

# ============================================================================
# Classes and Functions
# ============================================================================
class Nuclide:
    """" Class for loading and passing data """
    def __init__(self, N, A, M, B, D=True):
        self.N = N  # Nuclide Name, ex. H1, O16, U235, U238
        self.A = A  # Atomic Number
        self.M = M  # Atomic Mass
        self.B = B  # Tuple of Bools for [Elastic Scatter, Total, Fission]
                    # Ex. 2, U238 B = [1, 1, 1]
                    # Ex. 1, H1   B = [1, 1, 0], No fission data for Hydrogen.
        self.D = D  # This flag specifies whether the class is currenlty being
                    # used to create new data, i.e. needs BNL or Doppler data.
        self.s = {}  # Cross Section Data
        self.e = {}  # Energy Mesh Data  
        self.T = []  # Data Temperature Values

        print('='*79+'\nLoading BNL Data for '+self.N+'\n'+'='*79+'\n')
        if D:
            self.Load_BNL_Data()

    def Load_BNL_Data(self):
        """ Load plot Brookhaven National Laboratory data """        
        e, s = [[0], [0], [0]], [[0], [0], [0]]
        
        # Making filenames that will be searched for.
        Type = ['ES', 'Total', 'Fission']
        for i in range(3):
            if self.B[i]:
                name = 'Data/BNL/'+self.N+'_'+Type[i]+'.txt'
                e[i], s[i] = np.loadtxt(name, unpack=True, skiprows=1, 
                                        delimiter=',')

        # BNL files sometimes contain repeated energy-values. This causes zero
        # to appear in the denominator of the Doppler-Broadening equation. All
        # but the first cross-section and non-unique energy-value are removed.
        print('Removing non-unique energy mesh points and corresponding ' + \
              'cross sections...\n')
        for i in range(3):
            Indices = self.Unique_Indices(e[i])
            e[i] = np.array([e[i][j] for j in Indices])
            s[i] = np.array([s[i][j] for j in Indices])
        # BNL data is always evaluated at 300 deg. Kelvin.
        self.e[300], self.s[300] = e, s
        self.T.append(300)
        print('Finished saving unique BNL data to object...\n')

    def Unique_Indices(self, l):
        """" Return indices of unique and first cases of non-unique values"""
        seen = []
        res  = []
        for i, n in enumerate(l):
            if n not in seen:
                res.append(i)
                seen.append(n)
        return res

    def Load_Doppler_Data(self, Temperature):
        """ Load Doppler-Broadened cross sections """
        types = ['ES', 'Total', 'Fission']
        print('Loading Doppler-Broadened Data at '+str(Temperature)+
              ' for '+self.N+'...\n')
        for i in Temperature:
            e = [[0], [0], [0]]
            s = [[0], [0], [0]]
            for j in range(3):
                if self.B[j]:
                    name = 'Data/Doppler/' + self.N + '_'\
                         + types[j] + '_' + str(i) + '.txt'
                    e[j], s[j] = np.loadtxt(name, unpack=True, 
                                              delimiter=',')
            self.e[i], self.s[i] = e, s
            self.T.append(i)

def Background_Cross_Section(N, s):
    """ 
    Calculates background cross section for the case that H1, O16, U235, and
    U238 are the only nuclides that the core is composed of.  NOTE this is it-
    self an approximation.
    
    N = [N_H1, N_O16, N_U235, N_U238]
    s = [s_H1_es, s_O16_es, s_U235_potential, s_U238 not needed actually]
    """
    # Hmm... actually we only need three of the cross sections so will leave 
    # the last value in "s" as none for now.
    assert len(N) == len(s), \
    'There must be a number density "N" for each cross-section "s".'
    
    num = N[0]*s[0] + N[1]*s[1] + N[2]*s[2]
    den = N[3]
    return num/den


# ============================================================================
# Testing
# ============================================================================
if __name__ == '__main__':
    print('TESTING the project utilities classes and functions.\n')
    os.chdir('../')

    H1 = Nuclide('H1', 1, 1.008, [1, 1, 0])
    H1.Load_Doppler_Data([600, 900, 1200])
    
    O16 = Nuclide('O16', 16, 15.995, [1, 1, 0])
    O16.Load_Doppler_Data([600, 900, 1200])
    
    U235 = Nuclide('U235', 235, 235.044, [1, 1, 1])
    U235.Load_Doppler_Data([600, 900, 1200])
    
    U238 = Nuclide('U238', 238, 238.051, [1, 1, 1])
    U238.Load_Doppler_Data([600, 900, 1200])