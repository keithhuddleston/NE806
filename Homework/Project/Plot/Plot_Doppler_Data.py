""" Final Project for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
     
    Last Edit: Dec. 7, 2020
"""

# ============================================================================
# Import Statements
# ============================================================================
import matplotlib.pyplot as plt
import os

# Custom files written for project
os.chdir('../')
from Utilities.Utilities import Nuclide
            
# ============================================================================
# Function for plotting
# ============================================================================
def Plot_Data(N, T, Mode=1, Type=0):
    """ Plots data of the nuclide 'N' for given temperature 'T' """
    plt.rcParams.update({'font.size': 14})
    if Mode == 1:
        # Make sure data for specified temperature is loaded
        assert T in N.T, 'Data for Temperature "T" not loaded'
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
        ax.set_xscale('log'), ax.set_yscale('log')
        for i in range(len(N.B)):
            if N.B[i]:
                ax.plot(N.e[T][i], N.s[T][i])
        plt.xlabel('Energy [eV]'), plt.ylabel('Cross Section [Barns]')
        plt.title(N.N + ' Data at ' + str(T) + '[\N{DEGREE SIGN}K]')
        if N.B[-1] == 0:
            plt.legend(['Elastic Scatter', 'Total'])
        elif N.B[-1] == 1:
            plt.legend(['Elastic Scatter', 'Total', 'Fission'])
        plt.show()
    elif Mode == 2:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 3))
        ax.set_xscale('log'), ax.set_yscale('log')
        for i in N.T:
            ax.plot(N.e[i][Type], N.s[i][Type])
        plt.legend([r'$300\degree K$', r'$600\degree K$', 
                    r'$900\degree K$', r'$1200\degree K$'])
        plt.xlabel('Energy [eV]')
        plt.ylabel('Cross Section [Barns]')
        plt.grid(ls='--', axis='y')
        plt.xlim([6, 11])
        
# ============================================================================
# Load data and use plotting function
# ============================================================================
H1 = Nuclide('H1', 1, 1.008, [1, 1, 0])
H1.Load_Doppler_Data([600, 900, 1200])

O16  = Nuclide('O16', 16, 15.995, [1, 1, 0])
O16.Load_Doppler_Data([600, 900, 1200])

U235 = Nuclide('U235', 235, 235.044, [1, 1, 1])
U235.Load_Doppler_Data([600, 900, 1200])

U238 = Nuclide('U238', 238, 238.051, [1, 1, 1])
U238.Load_Doppler_Data([600, 900, 1200])

Plot_Data(H1,   1200, 2, 1)
Plot_Data(O16,  1200, 2, 0)
Plot_Data(U235, 1200, 2, 1)
Plot_Data(U238, 1200, 2, 1)