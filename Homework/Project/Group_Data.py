""" Final Project for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
     
    Create data files for the flux averaged cross sections. Note that plotting
    is being performed in this file as well, in contrast to the Doppler-Broad-
    ening step.
    
    Last Edit: Dec. 7, 2020
"""

# ============================================================================
# Import statements
# ============================================================================
import numpy as np

# Custom files written for project
from Utilities.Utilities import Nuclide


# ============================================================================
# Functions
# ============================================================================
def Group(e_vals, xs_vals, Group_Structure, e_total, xs_total, xs_dilution):
    """ Seperate input data into groups based on given group structure """
   
    e_vals, xs_vals = list(e_vals), list(xs_vals)
    e_total, xs_total = list(e_total), list(xs_total)
    
    nn = 1000000
    
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
    y_vals_group = []
    for i in range(shape):
        num = np.trapz(phi_group[i] * xs_group[i], e_group[i])
        den = np.trapz(phi_group[i], e_group[i])
        val = num/den
        y_vals_group.append(val)
        
    return y_vals_group[::-1]

def Make_Group_Data(Nuclide, Group_Structure, Dilution):
    """ Create flux weighted averages cross section data """
    Length = len(Group_Structure)
    Shape = Length - 1
    Folder_Name = 'Data/Group/'+str(Shape)+'/'
    print('='*79+'\nSaving '+Nuclide.N+' group data for '+\
          str(Shape)+' group structure to '+Folder_Name+'\n'+\
          '='*79+'\n')

    Temperatures = Nuclide.T
    data_type = ['ES', 'Total', 'Fission']
    for t in Temperatures:
        for i in range(len(Nuclide.B)):
            # Total cross-section data is always index 1
            et = Nuclide.e[t][1]
            st = Nuclide.s[t][1]
            if Nuclide.B[i]:
                print('Grouping data at '+str(t)+u"\N{DEGREE SIGN}K"+' for '
                      +data_type[i]+' XS\n')
                group_data = np.zeros((Length, len(Dilution)))
                for j in range(len(Dilution)):
                    group_data[0][j] = Dilution[j]
                    e = Nuclide.e[t][i]
                    s = Nuclide.s[t][i]
                    d = Dilution[j]
                    column = Group(e, s, Group_Structure, et, st, d)
                    for k in range(len(column)):
                        group_data[k+1][j] = column[k]
                File_Name = Folder_Name+Nuclide.N+'G_'+data_type[i]+'_'+str(t)+'.txt'
                np.savetxt(File_Name, group_data)
    # Nothing to return, all data is stored in text files.
    return


# ============================================================================
# Testing
# ============================================================================
if __name__ == '__main__':
    # --- Function Testing Section ---
    print('Testing the file Group_Data.py\n')
    
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
                         3.00e-8,  1.00e-11])*1e6# eV
    Casmo_16 = Casmo_16[::-1]
    
    Casmo_2 = np.array([1.00e1, 1.00e-7, 1.00e-11]) # MeV
    Casmo_2 = Casmo_2[::-1]

    # Microscopic Dilution/Background Cross Section
    Dilution = [1e1, 1e2, 1e3, 1e4, 1e5] # Barns

    # Test
    # e_vals, xs_vals, Group_Structure, e_total, xs_total, xs_dilution
    e_vals, xs_vals = H1.e[300][1], H1.s[300][1]
    e_total, xs_total = H1.e[300][1], H1.s[300][1]
    Test_0 = Group(e_vals, xs_vals, Casmo_16, e_total, xs_total, 101)
    
    Make_Group_Data(H1,   Casmo_2, Dilution)
    Make_Group_Data(O16,  Casmo_2, Dilution)
    Make_Group_Data(U235, Casmo_2, Dilution)
    Make_Group_Data(U238, Casmo_2, Dilution)
    
    Make_Group_Data(H1,   Casmo_16, Dilution)
    Make_Group_Data(O16,  Casmo_16, Dilution)
    Make_Group_Data(U235, Casmo_16, Dilution)
    Make_Group_Data(U238, Casmo_16, Dilution)