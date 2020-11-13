# ============================================================================
# Import statements
# ============================================================================
import matplotlib.pyplot as plt
import numpy as np
from NE806_Functions import Nuclide_Data # File written for class
from NE806_Functions import Seperate_Groups

# ============================================================================
# Functions
# ============================================================================
Es = np.array([1.00e1 , 8.21e-1, 5.53e-3, 4.00e-6, 1.30e-6, 1.15e-6, 1.097e-6,
               1.02e-6, 9.71e-7, 8.50e-7, 6.25e-7, 3.50e-7, 2.80e-7, 1.40e-7 ,
               5.80e-8, 3.00e-8, 1.00e-11])*1e6 # eV

def f(E, alpha, Ep):
    v = 1/Ep/(1-alpha)*np.ones_like(E)
    v[alpha*Ep >= E] = 0.0
    v[E >= Ep] = 0.0
    return v
R = np.zeros((len(Es)-1,len(Es)-1))
S = R*0
alpha = 0.0



nn = 1000
E_log = [np.linspace(Es[i], Es[i-1], nn) for i in range(1, len(Es))]
s = [np.zeros(nn)+i/len(Es)+1 for i in np.arange(len(Es))]
p = [1/i for i in E_log]

for g in range(len(Es)-1): # Minus 1 because Es are the bounds
    E_g = E_log[g]
    for gp in range(len(Es)-1):
        E_gp = E_log[gp]
        vals = [np.trapz(f(E_g, alpha, E_gp[i]), E_g) for i in range(len(E_gp))]
        R_gp_g = np.trapz(s[gp]*p[gp]*np.array(vals), E_gp)
        R[g, gp] = R_gp_g
        S[g, gp] = R[g, gp] / np.trapz(p[gp], E_gp)  

plt.matshow(S)
a = S[0]
for i in range(1, len(S)):
    a += S[i]