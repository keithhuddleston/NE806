""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
"""

# ============================================================================
# Import statements
# ============================================================================
import numpy as np
import matplotlib.pyplot as plt
from NE806_Functions import Nuclide_Data # File written for class


# ============================================================================
# Functions
# ============================================================================

# Narrow Resonance Flux Approximation.
def phinr(sigma_t, sigma_o, E):
    return 1 / (E*(sigma_t + sigma_o))

# Microscopic Flux Weighted Average Cross Section.  Note, that the energy mesh
# must be ultra fine for this simplifaction to work properly.
def sigma_g(Phi, XS):
    return(sum(Phi*XS)/sum(Phi))

# Function for seperating the cross section data obtained from either the BNL
# website or cross section data obtained using the Doppler function.  Note,
# 
def Seperate_Groups(x_vals, y_vals, group_structure):
    # Indices used to place group stucture energy bounds
    ind_G = np.zeros(len(group_structure), dtype=int)
    k = 0
    for i in range(len(x_vals)):
        if x_vals[i] > group_structure[k]:
            ind_G[k] = int(i)
            k += 1
    if ind_G[-1] == 0:
        ind_G[-1] = len(x_vals)

    xg = [] # Group-wise x values
    yg = [] # Group-wise y values
    
    # Interpolate group boundary values for y
    bg = np.interp(group_structure, x_vals, y_vals)
    
    # Seperate data into groupwise lists
    for i in range(len(ind_G[:-1])):
        L1 = [group_structure[i]] + list(x_vals[ind_G[i]:ind_G[i+1]]) \
           + [group_structure[i+1]]
        L2 = [bg[i]] + list(y_vals[ind_G[i]:ind_G[i+1]]) \
           + [bg[i+1]]
        if L1[0] == L1[1]:
            L1 = L1[1:]
            L2 = L2[1:]
        if L1[-1] == L1[-2]:
            L1 = L1[:-1]
            L2 = L2[:-1]            
        xg.append(np.array(L1))
        yg.append(np.array(L2))
    # Return a list containing lists of groupwise data
    return xg, yg

def Make_Group_Datafile(ic, sd, gs, filename, rt=1, Plot=True):
    Temps = [str(i) for i in ic.Temps]
    header = 'NG='+str(len(gs)-1)+'\nGL\tGU\tTemperature\t'
    for i in sd:
        header = header + str(i) + '\t'
    GL = gs[:-1]
    GU = gs[1:]
    N = np.zeros((len(Temps)*len(GU), len(sd)+3))
    for i in range(len(Temps)):
        for j in range(len(sd)):
            Emesh   = ic.ESEM[Temps[i]]
            sigma_a = ic.NGXS[Temps[i]]
            if rt == 1:
                xs = ic.ESXS[Temps[i]]
            elif rt ==2:
                xs = ic.NGXS[Temps[i]]
            sigma_e = ic.ESXS[Temps[i]]
            phi     = phinr(sigma_a+sigma_e, sd[j], Emesh)
            e, s = Seperate_Groups(Emesh, xs, gs)
            p = Seperate_Groups(Emesh, phi, gs)[1]
            sg = np.zeros(len(e))
            for k in range(len(sg)):
                sg[k] = sigma_g(p[k], s[k])
                N[k+i*len(sg)][2+j] = sg[k]
                N[k+i*len(sg)][2] = Temps[i]
                N[k+i*len(sg)][0]= GL[k]
                N[k+i*len(sg)][1]= GU[k]
            # Plot Results
            if Plot:
                fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
                ax1.set_xscale("log")
                ax1.set_yscale("log")
                plt.xlabel('Energy (eV)'), plt.ylabel('Cross Section (Barns)')
                ax1.plot(Emesh, xs, color='black', linewidth=2.25)
                plt.xlim(min(gs), max(gs)) 
                for k in range(len(gs)-1):
                    ax1.plot([gs[k], gs[k+1]], [sg[k], sg[k]], c='r', linewidth=2.25)
                for k in range(len(gs)-2):
                    ax1.plot([gs[k+1], gs[k+1]], [sg[k], sg[k+1]], c='r', linewidth=2.25)
                ax2 = ax1.twinx()
                ax2.loglog(Emesh, phi, linewidth=2.25)
                ax2.set_ylabel('Flux Shape Function Arb. Units')
                plt.legend([''])
                plt.show()
    np.savetxt(filename, N, delimiter='\t', header=header)
    return N


# ============================================================================
# Define Group and Background Data
# ============================================================================
# Casmo 2 group structure
C2 = np.array([1.00e1, 1.00e-7, 1.00e-11])*1e6 # eV
C2 = C2[::-1]
# Casmo 16 group structure
C16 = np.array([1.00e1 , 8.21e-1, 5.53e-3, 4.00e-6, 1.30e-6, 1.15e-6, 1.097e-6,
                1.02e-6, 9.71e-7, 8.50e-7, 6.25e-7, 3.50e-7, 2.80e-7, 1.40e-7 ,
                5.80e-8, 3.00e-8, 1.00e-11])*1e6 # eV
C16 = C16[::-1]

# Microscopic Dilution/Background Cross Section
sd = [1e1, 1e2, 1e3, 1e4, 1e5]


# ============================================================================
# Load Interpolated Interpreted Plotted Data Files and Doppler-Broadened Data
# ============================================================================
H1 = Nuclide_Data(1.008, False)
O16 = Nuclide_Data(15.995, False)
U238 = Nuclide_Data(238.051, True)
U235 = Nuclide_Data(235.044, True)

H1.load_data('Data/Doppler/H1_ES_300.txt', 'ES', 300)
H1.load_data('Data/Doppler/H1_ES_600',     'ES', 600)
H1.load_data('Data/Doppler/H1_ES_900',     'ES', 900)
H1.load_data('Data/Doppler/H1_ES_1200',    'ES', 1200)

O16.load_data('Data/Doppler/O16_ES_300.txt', 'ES', 300)
O16.load_data('Data/Doppler/O16_ES_600',     'ES', 600)
O16.load_data('Data/Doppler/O16_ES_900',     'ES', 900)
O16.load_data('Data/Doppler/O16_ES_1200',    'ES', 1200)

U238.load_data('Data/Doppler/U238_ES_300.txt', 'ES', 300)
U238.load_data('Data/Doppler/U238_ES_600',     'ES', 600)
U238.load_data('Data/Doppler/U238_ES_900',     'ES', 900)
U238.load_data('Data/Doppler/U238_ES_1200',    'ES', 1200)
U238.load_data('Data/Doppler/U238_NG_300.txt', 'NG', 300)
U238.load_data('Data/Doppler/U238_NG_600',     'NG', 600)
U238.load_data('Data/Doppler/U238_NG_900',     'NG', 900)
U238.load_data('Data/Doppler/U238_NG_1200',    'NG', 1200)

U235.load_data('Data/Doppler/U235_ES_300.txt', 'ES', 300)
U235.load_data('Data/Doppler/U235_ES_600',     'ES', 600)
U235.load_data('Data/Doppler/U235_ES_900',     'ES', 900)
U235.load_data('Data/Doppler/U235_ES_1200',    'ES', 1200)
U235.load_data('Data/Doppler/U235_NG_300.txt', 'NG', 300)
U235.load_data('Data/Doppler/U235_NG_600',     'NG', 600)
U235.load_data('Data/Doppler/U235_NG_900',     'NG', 900)
U235.load_data('Data/Doppler/U235_NG_1200',    'NG', 1200)


# ============================================================================
# Create Group Files
# ============================================================================
# H1GE   = Make_Group_Datafile(H1,   sd, C16, 'Data/Group/H1_Group_ES.txt',   1)
O16GE  = Make_Group_Datafile(O16,  sd, C16, 'Data/Group/O16_Group_ES.txt',  1)
# U238GE = Make_Group_Datafile(U238, sd, C16, 'Data/Group/U238_Group_ES.txt', 1)
# U238GN = Make_Group_Datafile(U238, sd, C16, 'Data/Group/U238_Group_NG.txt', 2)
# U235GE = Make_Group_Datafile(U235, sd, C16, 'Data/Group/U235_Group_ES.txt', 1)
# U235GB = Make_Group_Datafile(U235, sd, C16, 'Data/Group/U235_Group_NG.txt', 2)
