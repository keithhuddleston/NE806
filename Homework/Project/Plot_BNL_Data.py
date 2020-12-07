""" Final Project for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
     
    Plot data contained in the text files from Brookhaven National Lab website
    Example Link: https://www.nndc.bnl.gov/sigma/getPlotData.jsp
    
    Last Edit: Dec. 7, 2020
"""

# ============================================================================
# Import Statements
# ============================================================================
import matplotlib.pyplot as plt

# Custom files written for project
from Utilities.Utilities import Nuclide


# ============================================================================
# Function for plotting
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
            ax.plot(N.e[300][i], N.s[300][i])
    plt.xlabel('Energy [eV]'), plt.ylabel('Cross Section [Barns]')
    plt.title(N.N + ' Data at ' + str(T) + '[\N{DEGREE SIGN}K]')
    if N.B[-1] == 0:
        plt.legend(['Elastic Scatter', 'Total'])
    elif N.B[-1] == 1:
        plt.legend(['Elastic Scatter', 'Total', 'Fission'])
    plt.show()


# ============================================================================
# Load data and use plotting function
# ============================================================================
print('PLOTTING cross-section data obtained from the BNL website\n')
H1   = Nuclide('H1',   1,   1.008,   [1, 1, 0])
O16  = Nuclide('O16',  16,  15.995,  [1, 1, 0])
U235 = Nuclide('U235', 235, 235.044, [1, 1, 1])
U238 = Nuclide('U238', 238, 238.051, [1, 1, 1])

Plot_Data(H1,   300)
Plot_Data(O16,  300)
Plot_Data(U235, 300)
Plot_Data(U238, 300)