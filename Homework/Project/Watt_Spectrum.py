""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
"""

# ============================================================================
# Import statements
# ============================================================================
import numpy as np
import matplotlib.pyplot as plt


# ============================================================================
# Functions
# ============================================================================
def Watt_Spectrum(E):
    return 0.4865*np.sinh(np.sqrt(2*E))*np.exp(-E)

c16 = np.array([1.00e1 , 8.21e-1, 5.53e-3, 4.00e-6, 1.30e-6, 1.15e-6, 1.097e-6,
                1.02e-6, 9.71e-7, 8.50e-7, 6.25e-7, 3.50e-7, 2.80e-7, 1.40e-7 ,
                5.80e-8, 3.00e-8, 1.00e-11]) # MeV
c2 = np.array([1.00e1, 1.00e-7, 1.00e-11]) # MeV
c16 = c16[::-1]
c2 = c2[::-1]


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

# ============================================================================
# Testing
# ============================================================================
x = np.logspace(-5, 7, 75000)/10**6
y = Watt_Spectrum(x)

xg, yg = Seperate_Groups(x, y, c16)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
for i in range(len(xg)):
    ax.fill_between(xg[i], np.trapz(yg[i], xg[i]))
ax.semilogx(x,y,c='k',ls='--')
plt.xlabel('Fission Neutron Energy [MeV]')
plt.ylabel('Probability')

xg, yg = Seperate_Groups(x, y, c2)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 6))
for i in range(len(xg)):
    ax.fill_between(xg[i], np.trapz(yg[i], xg[i]))
ax.semilogx(x,y,c='k',ls='--')
plt.xlabel('Fission Neutron Energy [MeV]')
plt.ylabel('Probability')