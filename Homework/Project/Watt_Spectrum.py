""" Project for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
     
    This file contains function which will be used to provide the chi-spectrum
    matrix for groups following:
        
         1      2         n
    1 [chi_1, chi_1 ... chi_1]
    2 [chi_2, chi_2 ... chi_2]
    .
    .
    .
    n [chi_n, chi_n ... chi_n]
    
    this matrix is a numpy array.
    
    See Duderstadt and Hamilton, Nuclear Reactor Analysis, pg. 294 for the
    matrix interpretation of the eigenvalue problem.
"""

# ============================================================================
# Import statements
# ============================================================================
import matplotlib.pyplot as plt
import matplotlib.lines as Mlines
import numpy as np


# ============================================================================
# Functions
# ============================================================================
def Watt_Spectrum(E):
    """ Watt Spectrum as a function of energy in [MeV] """
    return 0.4865*np.sinh(np.sqrt(2*E))*np.exp(-E)

# Chi_Matrix is the function which will be imported to other Python files.
def Chi_Matrix(group_structure):
    """ Provides the Chi Matrix Component as described in Nuclear Reactor An-
        alysis 1976.
    """
    print('Computing the Chi-spectrum for a '+str(len(group_structure)-1)+\
          ' group structure...\n')
    shape = len(group_structure) - 1
    energy_group = []
    chi_group = []

    for i in range(1, shape+1):
        start = group_structure[i]
        end = group_structure[i-1]
        energy = np.linspace(start, end, 100)
        
        energy_group.append(energy)
        chi_group.append(Watt_Spectrum(energy))

    for i in range(shape):
        chi_group[i] = np.trapz(chi_group[i], energy_group[i])

    chi_matrix = np.zeros((shape, shape))

    for i in range(shape):
        chi_matrix[i] = np.zeros(shape)+chi_group[i]

    return chi_matrix


# ============================================================================
# Testing and Plotting
# ============================================================================
if __name__ == '__main__':
    # --- Function Testing Section ---
    Casmo_16 = np.array([1.00e1,  8.21e-1,  5.53e-3, 4.00e-6, 1.30e-6, 
                         1.15e-6, 1.097e-6, 1.02e-6, 9.71e-7, 8.50e-7, 
                         6.25e-7, 3.50e-7,  2.80e-7, 1.40e-7, 5.80e-8, 
                         3.00e-8, 1.00e-11]) # MeV

    Casmo_2 = np.array([1.00e1, 1.00e-7, 1.00e-11]) # MeV

    # Microscopic Dilution/Background Cross Section
    Dilution = [1e1, 1e2, 1e3, 1e4, 1e5] # Units [Barns]

    Chi_16 = Chi_Matrix(Casmo_16)
    Chi_2 = Chi_Matrix(Casmo_2)

    # --- Plotting Section ---
    plt.rcParams.update({'font.size': 15})
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 5))
    
    chi_group_16 = [i[0] for i in Chi_16]
    chi_group_2 = [i[0] for i in Chi_2]
    
    # Calculate Watt Spectrum
    e = np.logspace(-5, 7, 10000)/10**6
    chi = Watt_Spectrum(e)
    ax.loglog(e, chi, c='k')

    # Casmo 16 Group Plotting
    shape = len(Casmo_16) - 1
    for i in range(shape):
        ax.plot([Casmo_16[i], Casmo_16[i+1]], 
                [chi_group_16[i], chi_group_16[i]], c='r')
        for i in range(len(Casmo_16)-2):
                ax.plot([Casmo_16[i+1], Casmo_16[i+1]], 
                         [chi_group_16[i], chi_group_16[i+1]], c='r')
    
    # Casmo 2 Group Plotting
    shape = len(Casmo_2) - 1
    for i in range(shape):
        ax.plot([Casmo_2[i], Casmo_2[i+1]], 
                [chi_group_2[i], chi_group_2[i]], c='g')
        for i in range(len(Casmo_2)-2):
                ax.plot([Casmo_2[i+1], Casmo_2[i+1]], 
                         [chi_group_2[i], chi_group_2[i+1]], c='g')

    ax.set_ylabel('Probability, Arbitrary Units')
    ax.set_xlabel('Energy of Neutron Born from Fission [MeV]')
    
    # Legend, because there are many seperate lines we need to define handles
    line1 = Mlines.Line2D([], [], color='k', label='Watt Spectrum')
    line2 = Mlines.Line2D([], [], color='r', label='Casmo 16')
    line3 = Mlines.Line2D([], [], color='g', label='Casmo 2')
    plt.legend(handles=[line1, line2, line3])
    
    plt.grid(axis='y', ls='--')
    
    # We want to set our own ticks so that we can capture the 100% or 1e0 line
    plt.yticks([1e-11, 1e-9, 1e-7, 1e-5, 1e-3, 1e-1, 1e0])
    
    plt.xlim(1.00e-11, 1.00e1)

    plt.show()