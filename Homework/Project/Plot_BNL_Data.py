""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
"""

# ============================================================================
# Import Statements
# ============================================================================
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})

class Nuclide_Data:
    def __init__(self, N, M, B):
        self.N = N    # Name of Nuclide, ex. H1, O16, U235, U238 
        self.M = M    # Atomic Mass of Nuclide
        self.B = B    # Tuple of Bools for [Elastic Scatter, Total, Fission]
                      # ex 1. H1 B = [1, 1, 0], ex2, U238 B = [1, 1, 1]
        self.XS = {}  # Cross Section Data
        self.EM = {}  # Energy Mesh Data       
        self.T = []   # Temperature Values of data
        
        print('='*25+'\nLoading BNL Data for ' + self.N + '\n' + '='*25 + '\n')
        self.Load_BNL_Data()
        
    def uniqueIndexes(self, l):
        seen = []
        res = []
        for i, n in enumerate(l):
            if n not in seen:
                res.append(i)
                seen.append(n)
        return res
    
    def Load_BNL_Data(self, EM=None):
        # Make sure that BNL data is not already loaded in, Note, all BNL data
        # from the plotted data found on the BNL website is at 300 deg K
        assert 300 not in self.T, 'BNL Data Already Loaded'
        self.T.append(300)
        
        # Create the file names of the data to look for
        file_names = []
        data_type = ['ES', 'Total', 'Fission']
        for i in range(3):
            if self.B[i]:
                file_names.append('Data/BNL/' + self.N + '_' + data_type[i] + '.txt')
        
        # Load in the data from the file names list
        EM = [[0], [0], [0]]
        XS = [[0], [0], [0]]
        for i in range(len(file_names)):
            EM[i], XS[i] = np.loadtxt(file_names[i], unpack=True, skiprows=1,
                                      delimiter=',')
        
        print('Removing non-unique energy mesh values...\n')
        i1 = self.uniqueIndexes(EM[0])
        i2 = self.uniqueIndexes(EM[1])
        i3 = self.uniqueIndexes(EM[2])
        EM = [[EM[0][i] for i in i1],
             [EM[1][i] for i in i2],
             [EM[2][i] for i in i3]]
        XS = [[XS[0][i] for i in i1],
             [XS[1][i] for i in i2],
             [XS[2][i] for i in i3]]
        
        print('Saving BNL data to object...\n')
        self.EM[300] = EM
        self.XS[300] = XS
        
    def Plot_Data(self, T):
        # Make sure data for specified temperature is loaded
        assert T in self.T, 'Data for Temperature "T" not loaded'
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
        ax.set_xscale('log'), ax.set_yscale('log')
        for i in range(len(self.B)):
            if self.B[i]:
                ax.plot(self.EM[300][i], self.XS[300][i])
        plt.xlabel('Energy [eV]'), plt.ylabel('Cross Section [Barns]')
        plt.title(self.N + ' Data at ' + str(T) + '[\N{DEGREE SIGN}K]')
        if self.B[-1] == 0:
            plt.legend(['Elastic Scatter', 'Total'])
        elif self.B[-1] == 1:
            plt.legend(['Elastic Scatter', 'Total', 'Fission'])

# ============================================================================
# Load Data
# ============================================================================
if __name__ == "__main__":
    H1 = Nuclide_Data('H1', 1.008, [1, 1, 0])
    O16 = Nuclide_Data('O16', 15.995, [1, 1, 0])
    U235 = Nuclide_Data('U235', 235.044, [1, 1, 1])
    U238 = Nuclide_Data('U238', 238.051, [1, 1, 1])
    
    H1.Plot_Data(300)
    O16.Plot_Data(300)
    U235.Plot_Data(300)
    U238.Plot_Data(300)