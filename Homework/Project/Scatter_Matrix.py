# ============================================================================
# Import statements
# ============================================================================
import matplotlib.pyplot as plt
import numpy as np
from Project_Utilities import Nuclide_Data # File written for class


# ============================================================================
# Functions
# ============================================================================

def f(E, alpha, Ep):
    v = 1/Ep/(1-alpha)*np.ones_like(E)
    v[alpha*Ep >= E] = 0.0
    v[E >= Ep] = 0.0
    return v

def Seperate_Group(Nuclide, group_structure, xs_dilution, Temperature):
    assert Temperature in Nuclide.T, "Specified temperature not loaded."
    # First let's get the scattering energy and cross-section data
    e_vals  = list(Nuclide.EM[Temperature][0])
    xs_vals  = list(Nuclide.XS[Temperature][0])
    # Second let's get the total energy and cross-section data
    e_total  = list(Nuclide.EM[Temperature][1])
    xs_total = list(Nuclide.XS[Temperature][1])
    
    e_vals_boundary = list(group_structure)*2
    xs_vals_boundary = list(np.interp(e_vals_boundary, e_vals, xs_vals))

    # To accound for the additional peaks introduced from the flux term the
    # points at the energy values of the total cross section are added.
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
            if bound_index == len(group_structure):
                break

    xs_total = np.interp(e_vals, e_total, xs_total)
    # Narrow Resonance Flux Approximation, see (ADD REFERENCE)
    phi = 1 / (e_vals * (xs_total + xs_dilution))

    e_group = []
    xs_group = []
    phi_group = []
    
    for i in range(len(indices)-1):
        i1 = indices[i]+1
        i2 = indices[i+1]+1
        e_group.append(e_vals[i1:i2])
        xs_group.append(xs_vals[i1:i2])
        phi_group.append(phi[i1:i2])
 
    e_group = e_group[::-1]
    xs_group = xs_group[::-1]
    phi_group = phi_group[::-1]
        
    return e_group, xs_group, phi_group

def Micro_Scatter_Matrix(Nuclide, group_structure, Dilution, Temperature):
    e, s, p = Seperate_Group(Nuclide, group_structure, Dilution, Temperature)
    alpha = ((Nuclide.A-1)/(Nuclide.A+1))**2
    R = np.zeros((len(group_structure)-1,len(group_structure)-1))
    S = R*0
    Es = group_structure[::-1]
    for g in range(len(Es)-1):
        E_g = e[g]
        for gp in range(len(Es)-1):
            E_gp = e[gp]
            vals = []
            for i in range(len(E_gp)):
                vals.append(np.trapz(f(E_g, alpha, E_gp[i]), E_g))
            R_gp_g = np.trapz(s[gp]*p[gp]*np.array(vals), E_gp)
            R[g, gp] = R_gp_g
            S[g, gp] = R[g, gp] / np.trapz(p[gp], E_gp)
    # plt.matshow(S)
    # a = S[0]
    # for i in range(1, len(S)):
    #     a += S[i]
    # print(a)
    return S
    
def Scatter_Matrix(Nuclides, Ns, group_structure, Dilution, Temperature):
    S = np.zeros((len(group_structure)-1, len(group_structure)-1))  
    for i in range(len(Nuclides)):
        S += S + Micro_Scatter_Matrix(Nuclides[i], group_structure, 
                                      Dilution, Temperature)*Ns[i]
    plt.matshow(S)
    return S

# ============================================================================
# Load Interpolated Interpreted Plotted Data Files and Doppler-Broadened Data
# ============================================================================
H1 = Nuclide_Data('H1', 1.008, [1, 1, 0], 1)
O16 = Nuclide_Data('O16', 15.995, [1, 1, 0], 16)
U235 = Nuclide_Data('U235', 235.044, [1, 1, 1], 235)
U238 = Nuclide_Data('U238', 238.051, [1, 1, 1], 238)

H1.Load_Doppler_Data([600, 900, 1200])
O16.Load_Doppler_Data([600, 900, 1200])
U235.Load_Doppler_Data([600, 900, 1200])
U238.Load_Doppler_Data([600, 900, 1200])

Casmo_16 = np.array([1.00e1,   8.21e-1, 5.53e-3, 4.00e-6, 1.30e-6, 1.15e-6, 
                     1.097e-6, 1.02e-6, 9.71e-7, 8.50e-7, 6.25e-7, 3.50e-7, 
                     2.80e-7, 1.40e-7,  5.80e-8, 3.00e-8, 1.00e-11])*1e6 # eV
Casmo_16 = Casmo_16[::-1]

# Microscopic Dilution/Background Cross Section
Dilution = [1e1, 1e2, 1e3, 1e4, 1e5] # Barns
# ============================================================================
# Testing
# ============================================================================
if __name__ == '__main__':
    Micro_Scatter_Matrix(H1, Casmo_16, Dilution[0], 300)
    Nuclides = [H1, O16, U235, U238]
    Ns = [6, 7, 0.1, 1.9]
    Scatter_Matrix(Nuclides, Ns, Casmo_16, Dilution[0], 300)