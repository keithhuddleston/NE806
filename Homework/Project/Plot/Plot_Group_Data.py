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
import matplotlib.pyplot as plt
import matplotlib.lines as Mlines
import numpy as np
import os


# Custom files written for project
os.chdir('../')
from Utilities.Utilities import Nuclide

def Group_Plot(Nuclide, Temperature, Type, group_structure):
    file_name = 'Data/Group/16/'+Nuclide.N+'G_'+Type+'_'+str(300)+'.txt'
    group = np.loadtxt(file_name, skiprows=1, unpack=True)    
    file_name = 'Data/BNL/'+Nuclide.N+'_'+Type+'.txt'
    e, s = np.loadtxt(file_name, skiprows=1, unpack=True, delimiter=',')
    file_name = 'Data/BNL/'+Nuclide.N+'_Total.txt'
    et, st = np.loadtxt(file_name, skiprows=1, unpack=True, delimiter=',')
    plt.rcParams.update({'font.size': 14})
    fig, ax1 = plt.subplots(nrows=2, ncols=1, figsize=(6, 8))
    ax1[0].grid(ls='--', axis='y')
    ax1[0].set_xscale("log")
    ax1[0].set_yscale("log")
    plt.xlabel('Energy (eV)'), ax1[0].set_ylabel('Cross Section (Barns)')
    ax1[0].plot(e, s, color='black')
    y_vals_group1 = group[0]
    y_vals_group2 = group[3]
    ax1[0].set_xlim(min(group_structure), max(group_structure)) 
    for k in range(len(group_structure)-1):
        ax1[0].plot([group_structure[k], group_structure[k+1]], 
                  [y_vals_group1[k], y_vals_group1[k]], c='r')
        for k in range(len(group_structure)-2):
            ax1[0].plot([group_structure[k+1], group_structure[k+1]], 
                      [y_vals_group1[k], y_vals_group1[k+1]], c='r')
    for k in range(len(group_structure)-1):
        ax1[0].plot([group_structure[k], group_structure[k+1]], 
                  [y_vals_group2[k], y_vals_group2[k]], c='g')
        for k in range(len(group_structure)-2):
            ax1[0].plot([group_structure[k+1], group_structure[k+1]], 
                      [y_vals_group2[k], y_vals_group2[k+1]], c='g')

    phi1 = 1 / (et * (st + 1))
    phi2 = 1 / (et * (st + 1000))

    line1 = Mlines.Line2D([], [], color='k', label=r'$\sigma_{t}(E)$')
    line2 = Mlines.Line2D([], [], color='r', label=r'$\sigma_{tg},\ d=1$')
    line3 = Mlines.Line2D([], [], color='g', label=r'$\sigma_{tg},\ d=10000$')
    line4 = Mlines.Line2D([], [], color='r',ls='--', label=r'$\phi (E),\ d=1$')
    line5 = Mlines.Line2D([], [], color='g',ls='--', label=r'$\phi (E),\ d=10000$')
    ax1[0].legend(handles=[line1, line2, line3])
    ax1[1].set_ylabel('Flux Shape')
    ax1[1].loglog(et, phi1, c='r', ls='--')
    ax1[1].loglog(et, phi2, c='g', ls='--')
    ax1[1].set_yticks([])
    ax1[1].set_xlim(min(group_structure), max(group_structure)) 
    ax1[0].set_xticks([])
    plt.subplots_adjust(hspace=.0)
    ax1[1].legend(handles=[line4, line5])
    return

def Temperature_Plot(Nuclide, Temperature, Type, group_structure):
    file_name = 'Data/Group/16/'+Nuclide.N+'G_'+Type+'_'+str(300)+'.txt'
    group = np.loadtxt(file_name, skiprows=1, unpack=True)    
    file_name = 'Data/BNL/'+Nuclide.N+'_'+Type+'.txt'
    e, s = np.loadtxt(file_name, skiprows=1, unpack=True, delimiter=',')
    file_name = 'Data/BNL/'+Nuclide.N+'_Total.txt'
    et, st = np.loadtxt(file_name, skiprows=1, unpack=True, delimiter=',')
    plt.rcParams.update({'font.size': 14})
    fig, ax1 = plt.subplots(nrows=2, ncols=1, figsize=(6, 8))
    ax1[0].grid(ls='--', axis='y')
    ax1[0].set_xscale("log")
    ax1[0].set_yscale("log")
    plt.xlabel('Energy (eV)'), ax1[0].set_ylabel('Cross Section (Barns)')
    ax1[0].plot(e, s, color='black')
    y_vals_group1 = group[0]
    y_vals_group2 = group[3]
    ax1[0].set_xlim(min(group_structure), max(group_structure)) 
    for k in range(len(group_structure)-1):
        ax1[0].plot([group_structure[k], group_structure[k+1]], 
                  [y_vals_group1[k], y_vals_group1[k]], c='r')
        for k in range(len(group_structure)-2):
            ax1[0].plot([group_structure[k+1], group_structure[k+1]], 
                      [y_vals_group1[k], y_vals_group1[k+1]], c='r')
    for k in range(len(group_structure)-1):
        ax1[0].plot([group_structure[k], group_structure[k+1]], 
                  [y_vals_group2[k], y_vals_group2[k]], c='g')
        for k in range(len(group_structure)-2):
            ax1[0].plot([group_structure[k+1], group_structure[k+1]], 
                      [y_vals_group2[k], y_vals_group2[k+1]], c='g')

    phi1 = 1 / (et * (st + 1))
    phi2 = 1 / (et * (st + 1000))

    line1 = Mlines.Line2D([], [], color='k', label=r'$\sigma_{t}(E)$')
    line2 = Mlines.Line2D([], [], color='r', label=r'$\sigma_{tg},\ d=1$')
    line3 = Mlines.Line2D([], [], color='g', label=r'$\sigma_{tg},\ d=10000$')
    line4 = Mlines.Line2D([], [], color='r',ls='--', label=r'$\phi (E),\ d=1$')
    line5 = Mlines.Line2D([], [], color='g',ls='--', label=r'$\phi (E),\ d=10000$')
    return

# ============================================================================
# Loading data, and defining parameters
# ============================================================================   
H1 = Nuclide('H1', 1, 1.008, [1, 1, 0])
H1.Load_Doppler_Data([600, 900, 1200])

# O16  = Nuclide('O16', 16, 15.995, [1, 1, 0])
# O16.Load_Doppler_Data([600, 900, 1200])

U235 = Nuclide('U235', 235, 235.044, [1, 1, 1])
U235.Load_Doppler_Data([600, 900, 1200])

# U238 = Nuclide('U238', 238, 238.051, [1, 1, 1])
# U238.Load_Doppler_Data([600, 900, 1200])

Casmo_16 = np.array([1.00e1,   8.21e-1,  5.53e-3, 4.00e-6, 1.30e-6, 
                     1.15e-6,  1.097e-6, 1.02e-6, 9.71e-7, 8.50e-7, 
                     6.25e-7,  3.50e-7,  2.80e-7, 1.40e-7, 5.80e-8, 
                     3.00e-8,  1.00e-11])*1e6 # eV
Casmo_2 = np.array([1.00e1, 1.00e-7, 1.00e-11])*1e6 # eV

# ============================================================================
# Plotting data
# ============================================================================
# Group_Plot(U235, 300, 'Total', Casmo_16)
Temperature_Plot(U235, 300, 'Total', Casmo_16)