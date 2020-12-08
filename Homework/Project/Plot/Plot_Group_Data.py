""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
     
    Create data files for the flux averaged cross sections. Note that plotting
    is being performed in this file as well, in contrast to the Doppler-Broad-
    ening step.
"""

# ============================================================================
# Import statements
# ============================================================================
import numpy as np
import matplotlib.pyplot as plt

# File written for class
from Project_Utilities import Nuclide_Data

def Group_Plot(Nuclide, Temperature, Type, Group_Structure):
    file_name = 'Data/Group/'+Nuclide.N+'G_'+Type+'_'+str(Temperature)+'.txt'
    data = np.loadtxt(file_name, skiprows=1, unpack=True)    

    # Nothing to Return
    return

    # plt.rcParams.update({'font.size': 14})
    # fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
    # ax1.set_xscale("log")
    # ax1.set_yscale("log")
    # plt.xlabel('Energy (eV)'), plt.ylabel('Cross Section (Barns)')
    # ax1.plot(e_vals, xs_vals, color='black')
    # plt.xlim(min(group_structure), max(group_structure)) 
    # for k in range(len(group_structure)-1):
    #     ax1.plot([group_structure[k], group_structure[k+1]], 
    #               [y_vals_group[k], y_vals_group[k]], c='r')
    #     for k in range(len(group_structure)-2):
    #         ax1.plot([group_structure[k+1], group_structure[k+1]], 
    #                   [y_vals_group[k], y_vals_group[k+1]], c='r')
    # ax2 = ax1.twinx()
    # ax2.loglog(e_vals, phi)
    # ax2.set_ylabel('Flux Shape Function Arb. Units')
    # plt.title(Name+' '+Type+' Cross-Section Evaluated at '
    #           +str(Temp)+r'$\degree\ K\ for\ $'+r'$\sigma_b=$'
    #           +str(xs_dilution)+' Barnes')
    
# ============================================================================
# Loading data, and defining parameters
# ============================================================================   
H1 = Nuclide_Data('H1', 1.008, [1, 1, 0], 1)
H1.Load_Doppler_Data([600, 900, 1200])

# O16 = Nuclide_Data('O16', 15.995, [1, 1, 0], 16)
# O16.Load_Doppler_Data([600, 900, 1200])

# U235 = Nuclide_Data('U235', 235.044, [1, 1, 1], 235)
# U235.Load_Doppler_Data([600, 900, 1200])

# U238 = Nuclide_Data('U238', 238.051, [1, 1, 1], 238)
# U238.Load_Doppler_Data([600, 900, 1200])

Casmo_16 = np.array([1.00e1,   8.21e-1,  5.53e-3, 4.00e-6, 1.30e-6, 
                     1.15e-6,  1.097e-6, 1.02e-6, 9.71e-7, 8.50e-7, 
                     6.25e-7,  3.50e-7,  2.80e-7, 1.40e-7, 5.80e-8, 
                     3.00e-8,  1.00e-11])*1e6 # eV

# ============================================================================
# Plotting data
# ============================================================================
Group_Plot(H1, 300, 'Total')