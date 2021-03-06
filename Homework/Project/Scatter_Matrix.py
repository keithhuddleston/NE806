# ============================================================================
# Import statements
# ============================================================================
import matplotlib.pyplot as plt
import numpy as np

# Files written for project
from Utilities.Utilities import Nuclide
from Utilities.Utilities import Background_Cross_Section

# ============================================================================
# Functions
# ============================================================================
def f(E, alpha, Ep):
    v = 1/Ep/(1-alpha)*np.ones_like(E)
    v[alpha*Ep >= E] = 0.0
    v[E >= Ep] = 0.0
    return v
    
def Seperate_Group(Nuclide, Group_Structure, xs_dilution, Temperature):
    assert Temperature in Nuclide.T, "Specified temperature not loaded."
    
    # First let's get the scattering energy and cross-section data
    e_vals  = list(Nuclide.e[Temperature][0])
    xs_vals  = list(Nuclide.s[Temperature][0])
    # Second let's get the total energy and cross-section data
    e_total  = list(Nuclide.e[Temperature][1])
    xs_total = list(Nuclide.s[Temperature][1])
    
    nn = 10000
    
    Length = len(Group_Structure)
    Indices = np.zeros(Length, dtype=int)
    for j in range(Length):
        for i in e_total[Indices[j-1]:]:
            if i <= Group_Structure[j]:
                Indices[j] = Indices[j] + 1
            else:
                break
        Indices[j] = Indices[j] + Indices[j-1]

    shape = len(Group_Structure)-1
    e_group = [0] * shape
    for i in range(shape):
        start = Group_Structure[i]
        end = Group_Structure[i+1]
        add = e_total[Indices[i]:Indices[i+1]]
        e_group[i] = list(np.linspace(start, end, nn)) + list(add)
    e_group = [np.sort(i) for i in e_group]
    
    xs_group = [np.interp(i, e_vals, xs_vals) for i in e_group]

    # Narrow Resonance Flux Approximation, see (ADD REFERENCE)
    # phi = 1 / (e_vals * (xs_total + xs_dilution))
    xs_total_group = [np.interp(i, e_total, xs_total) for i in e_group]
    phi_group = [1 / (e_group[i] * (xs_total_group[i] + xs_dilution)) for i \
                  in range(shape)]
    
    e_group = e_group[::-1]
    xs_group = xs_group[::-1]
    phi_group = phi_group[::-1]
        
    return e_group, xs_group, phi_group

def Micro_Scatter_Matrix(Nuclide, group_structure, Dilution, Temperature):
    e, s, p = Seperate_Group(Nuclide, group_structure, Dilution, Temperature)
    alpha = ((Nuclide.A-1)/(Nuclide.A+1))**2
    shape = len(group_structure)-1
    R = np.zeros((shape, shape))
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
        matrix = Micro_Scatter_Matrix(Nuclides[i], group_structure, 
                                      Dilution, Temperature)
        S = S + matrix*Ns[i]
    plt.matshow(S)
    return S

# ============================================================================
# Testing
# ============================================================================
if __name__ == '__main__':
    print('Testing the file Scatter_Matrix.py\n')
    # Function Input
    N = [6, 7, 0.1, 1.9]
    s = [20, 10, 10, 0]
    
    H1 = Nuclide('H1', 1, 1.008, [1, 1, 0])
    H1.Load_Doppler_Data([600, 900, 1200])
    
    O16 = Nuclide('O16', 16, 15.995, [1, 1, 0])
    O16.Load_Doppler_Data([600, 900, 1200])
    
    U235 = Nuclide('U235', 235, 235.044, [1, 1, 1])
    U235.Load_Doppler_Data([600, 900, 1200])
    
    U238 = Nuclide('U238', 238, 238.051, [1, 1, 1])
    U238.Load_Doppler_Data([600, 900, 1200])

    Casmo_16 = np.array([1.00e1,   8.21e-1,  5.53e-3, 4.00e-6, 1.30e-6, 
                         1.15e-6,  1.097e-6, 1.02e-6, 9.71e-7, 8.50e-7, 
                         6.25e-7,  3.50e-7,  2.80e-7, 1.40e-7, 5.80e-8, 
                         3.00e-8,  1.00e-11])*1e6 # eV
    Casmo_16 = Casmo_16[::-1]
    
    Casmo_2 = np.array([1.00e1, 1.00e-7, 1.00e-11])*1e6 # eV

    so = Background_Cross_Section(N, s)

    Seperate_Group(H1, Casmo_16, so, 300)

    Test_1 = Micro_Scatter_Matrix(O16, Casmo_16, 100, 300)
    a = Test_1[0]
    for i in range(1, len(a)):
        a += Test_1[i]
    print(a)
    
    Nuclides = [H1, O16, U235, U238]
    Test_3 = Scatter_Matrix(Nuclides, N, Casmo_16, so, 300)
    