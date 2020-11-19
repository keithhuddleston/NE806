""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
"""

# ============================================================================
# Import Statements
# ============================================================================
import numpy as np

class Nuclide_Data:
    """" Class for loading and passing data """
    def __init__(self, N, M, B):
        self.N = N    # Name of Nuclide, ex. H1, O16, U235, U238 
        self.M = M    # Atomic Mass of Nuclide
        self.B = B    # Tuple of Bools for [Elastic Scatter, Total, Fission]
                      # ex 1. H1 B = [1, 1, 0], ex2, U238 B = [1, 1, 1]
        self.XS = {}  # Cross Section Data
        self.EM = {}  # Energy Mesh Data       
        self.T = []   # Temperature Values of data
        
        print('='*25+'\nLoading BNL Data for '+self.N+'\n'+'='*25+'\n')
        self.Load_BNL_Data()

    def Load_BNL_Data(self, EM=None):
        """ Load plot data from Brookhaven National Laboratory """
        self.T.append(300)
        
        EM, XS = [[0], [0], [0]], [[0], [0], [0]]
        types = ['ES', 'Total', 'Fission']
        for i in range(3):
            if self.B[i]:
                name = 'Data/BNL/'+self.N+'_'+types[i]+'.txt'
                EM[i], XS[i] = np.loadtxt(name, unpack=True, skiprows=1, 
                                          delimiter=',')
    
        # BNL files sometimes contain repeated energy-values. This causes zero
        # to appear in the denominator of the Doppler-Broadening equation. All
        # but the first cross-section and non-unique energy-value are removed.
        print('Removing non-unique EM and corresponding XS values...\n')
        i1 = self.uniqueIndexes(EM[0])
        i2 = self.uniqueIndexes(EM[1])
        i3 = self.uniqueIndexes(EM[2])
        EM = [np.array([EM[0][i] for i in i1]),
             np.array([EM[1][i] for i in i2]),
             np.array([EM[2][i] for i in i3])]
        XS = [np.array([XS[0][i] for i in i1]),
             np.array([XS[1][i] for i in i2]),
             np.array([XS[2][i] for i in i3])]
        
        print('Saving BNL data to object...\n')
        self.EM[300], self.XS[300] = EM, XS

    def uniqueIndexes(self, l):
        """" Return indices of unique and first cases of non-unique values"""
        seen = []
        res = []
        for i, n in enumerate(l):
            if n not in seen:
                res.append(i)
                seen.append(n)
        return res

    def Load_Doppler_Data(self, Temp):
        """ Load Doppler-Broadened cross sections """
        types = ['ES', 'Total', 'Fission']
        print('Loading Doppler-Broadened Data...')
        for i in Temp:
            EM = [[0], [0], [0]]
            XS = [[0], [0], [0]]
            for j in range(3):
                if self.B[j]:
                    name = 'Data/Doppler/' + self.N + '_'\
                         + types[j] + '_' + str(i) + '.txt'
                    EM[j], XS[j] = np.loadtxt(name, unpack=True, delimiter=',')
            self.EM[i], self.XS[i] = EM, XS
            self.T.append(i)