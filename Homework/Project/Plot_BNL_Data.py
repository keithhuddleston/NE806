""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
     
    Plot data contained in the text files from Brookhaven National Lab website
    Example Link: https://www.nndc.bnl.gov/sigma/getPlotData.jsp
"""

# ============================================================================
# Import Statements
# ============================================================================
import matplotlib.pyplot as plt

from Project_Utilities import Nuclide_Data

# ============================================================================
# Load and Plot Data
# ============================================================================
def Plot_Data(N, T):
    """ Plots data of the nuclide 'N' for given temperature 'T' """
    plt.rcParams.update({'font.size': 18})
    # Make sure data for specified temperature is loaded
    assert T in N.T, 'Data for Temperature "T" not loaded'
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
    ax.set_xscale('log'), ax.set_yscale('log')
    for i in range(len(N.B)):
        if N.B[i]:
            ax.plot(N.EM[300][i], N.XS[300][i])
    plt.xlabel('Energy [eV]'), plt.ylabel('Cross Section [Barns]')
    plt.title(N.N + ' Data at ' + str(T) + '[\N{DEGREE SIGN}K]')
    if N.B[-1] == 0:
        plt.legend(['Elastic Scatter', 'Total'])
    elif N.B[-1] == 1:
        plt.legend(['Elastic Scatter', 'Total', 'Fission'])
    plt.show()

if __name__ == "__main__":
    H1 = Nuclide_Data('H1', 1.008, [1, 1, 0])
    O16 = Nuclide_Data('O16', 15.995, [1, 1, 0])
    U235 = Nuclide_Data('U235', 235.044, [1, 1, 1])
    U238 = Nuclide_Data('U238', 238.051, [1, 1, 1])
    
    Plot_Data(H1, 300)
    Plot_Data(O16, 300)
    Plot_Data(U235, 300)
    Plot_Data(U238, 300)