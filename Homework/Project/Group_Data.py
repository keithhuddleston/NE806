""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
     
    Create data files for the flux averaged cross sections
"""

# ============================================================================
# Import statements
# ============================================================================
import numpy as np
import matplotlib.pyplot as plt
from Project_Utilities import Nuclide_Data # File written for class


# ============================================================================
# Functions
# ============================================================================
def Group(e_vals, xs_vals, group_structure, e_total, xs_total, xs_dilution, total=False, plot=False):
    """ Seperate input data into groups based on given group structure """
    print("="*34+'\nCalculating Group Cross Section...\n'+'='*34+'\n')
   
    e_vals, xs_vals = list(e_vals), list(xs_vals)
    e_total, xs_total = list(e_total), list(xs_total)
    e_vals_boundary = list(group_structure)*2
    xs_vals_boundary = list(np.interp(e_vals_boundary, e_vals, xs_vals))
    # When grouping the total cross section only boundary values need be added
    if total:
        xs_vals = xs_vals + xs_vals_boundary
        e_vals = e_vals + e_vals_boundary
    # To accound for the additional peaks introduced from the flux term the
    # points at the energy values of the total cross section are added.
    else:
        xs_vals_new = list(np.interp(e_total, e_vals, xs_vals))
        xs_vals = xs_vals + xs_vals_boundary + xs_vals_new
        e_vals = e_vals + e_vals_boundary + e_total
    exs_sorted = sorted(zip(e_vals, xs_vals))
    xs_vals = [xs for e, xs in exs_sorted]
    e_vals = [e for e, xs in exs_sorted]
    
    indices = np.zeros(len(group_structure), dtype=int)
    bound_index = 0
    for i, val in enumerate(e_vals):
        if val == group_structure[bound_index]:
            indices[bound_index] = i
            bound_index = bound_index + 1
            if bound_index == 17:
                break
            
    xs_total = np.interp(e_vals, e_total, xs_total)
    # Narrow Resonance Flux Approximation, see (ADD REFERENCE)
    phi = 1 / (e_vals * (xs_total + xs_dilution))

    y_vals_group = []
    for i in range(len(indices)-1):
        num = np.trapz(phi[indices[i]+1:indices[i+1]+1] \
                       * xs_vals[indices[i]+1:indices[i+1]+1], \
                       e_vals[indices[i]+1:indices[i+1]+1])
        den = np.trapz(phi[indices[i]+1:indices[i+1]+1], \
                       e_vals[indices[i]+1:indices[i+1]+1])
        val = num/den
        y_vals_group.append(val)
        
    if plot:
        plt.rcParams.update({'font.size': 14})
        fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        plt.xlabel('Energy (eV)'), plt.ylabel('Cross Section (Barns)')
        ax1.plot(e_vals, xs_vals, color='black')
        plt.xlim(min(group_structure), max(group_structure)) 
        for k in range(len(group_structure)-1):
            ax1.plot([group_structure[k], group_structure[k+1]], [y_vals_group[k], y_vals_group[k]], c='r')
            for k in range(len(group_structure)-2):
                ax1.plot([group_structure[k+1], group_structure[k+1]], [y_vals_group[k], y_vals_group[k+1]], c='r')
        ax2 = ax1.twinx()
        ax2.loglog(e_vals, phi)
        ax2.set_ylabel('Flux Shape Function Arb. Units')
    return y_vals_group

# def Make_Group_Data(N, Group_Structure, Dilution):
#     """ Create flux weighted averages cross section data """
#     Temperatures = N.T
#     for t in Temperatures:
#         for i in N.B:
#             if i:
#                 # Index '1' should always correspond to total cross section
#                 phi = Narrow_Resonance_Approx(N.XS[t][1], Dilution[0], 
#                                               N.EM[t][1], N.EM[t][i])
#                 phi = Seperate(N.EM[t][i], phi, Group_Structure)
#                 xs = Seperate(N.EM[t][i], N.XS[t][i])
#             xs_group = np.zeros(len(Group_Structure))
#             for j in range(len(Group_Structure)):
#                 xs_group[j] = sigma_g(phi[j], xs[j])
#     return


# ============================================================================
# Define Group and Background Data
# ============================================================================
Casmo_2 = np.array([1.00e1, 1.00e-7, 1.00e-11])*1e6 # eV
Casmo_2 = Casmo_2[::-1]

Casmo_16 = np.array([1.00e1,   8.21e-1, 5.53e-3, 4.00e-6, 1.30e-6, 1.15e-6, 
                     1.097e-6, 1.02e-6, 9.71e-7, 8.50e-7, 6.25e-7, 3.50e-7, 
                     2.80e-7, 1.40e-7,  5.80e-8, 3.00e-8, 1.00e-11])*1e6 # eV
Casmo_16 = Casmo_16[::-1]

# Microscopic Dilution/Background Cross Section
Dilution = [1e1, 1e2, 1e3, 1e4, 1e5]


# ============================================================================
# Load Interpolated Interpreted Plotted Data Files and Doppler-Broadened Data
# ============================================================================
H1 = Nuclide_Data('H1', 1.008, [1, 1, 0])
O16 = Nuclide_Data('O16', 15.995, [1, 1, 0])
U235 = Nuclide_Data('U235', 235.044, [1, 1, 1])
# U238 = Nuclide_Data('U238', 238.051, [1, 1, 1])

H1.Load_Doppler_Data([600, 900, 1200])
# O16.Load_Doppler_Data([600, 900, 1200])
# U235.Load_Doppler_Data([600, 900, 1200])
# U238.Load_Doppler_Data([600, 900, 1200])

x = U235.EM[300][0]
y = U235.XS[300][0]
xt = U235.EM[300][1]
yt = U235.XS[300][1]
d = Dilution[0]
Group(x, y, Casmo_16, xt, yt, d, False, True)

x = H1.EM[300][0]
y = H1.XS[300][0]
xt = H1.EM[300][1]
yt = H1.XS[300][1]
d = Dilution[0]
A = Group(x, y, Casmo_16, xt, yt, d, False, True)
# ============================================================================
# Create Group Files
# ============================================================================
# Make_Group_Data(H1, Casmo_16, Dilution)
