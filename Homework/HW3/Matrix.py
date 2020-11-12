# ============================================================================
# Import statements
# ============================================================================
import matplotlib.pyplot as plt
import numpy as np
from NE806_Functions import Nuclide_Data # File written for class
from NE806_Functions import Seperate_Groups
from NE806_Functions import phinr

# ============================================================================
# Load Interpolated Interpreted Plotted Data Files and Doppler-Broadened Data
# ============================================================================
H1 = Nuclide_Data(1.008, False)
H1.load_data('Data/Doppler/H1_ES_300.txt', 'ES', 300)
plt.loglog(H1.ESEM['300'], H1.ESXS['300'])

# ============================================================================
# Functions
# ============================================================================
Es = np.array([1.00e1 , 8.21e-1, 5.53e-3, 4.00e-6, 1.30e-6, 1.15e-6, 1.097e-6,
               1.02e-6, 9.71e-7, 8.50e-7, 6.25e-7, 3.50e-7, 2.80e-7, 1.40e-7 ,
               5.80e-8, 3.00e-8, 1.00e-11])*1e6 # eV
def Scatter_Matrix(sigma, N, A, phi):
    Sigma = sigma*N
    alpha = ((A-1)/(A+1))**2
def f(E, alpha, Ep):
    v = 1/Ep/(1-alpha)*np.ones_like(E)
    v[alpha*Ep >= E] = 0.0
    v[E >= Ep] = 0.0
    return v
R = np.zeros((len(Es)-1,len(Es)-1))
S = R*0

alpha = 0.0
# ============================================================================
# Functions
# ============================================================================
Em = H1.ESEM['300']
sigma_e = H1.ESXS['300']
sigma_a = H1.NGXS['300']
sd = [1e1, 1e2, 1e3, 1e4, 1e5]
phi_vals = phinr(sigma_a+sigma_e, sd[0], Em)
e, s = Seperate_Groups(Em, sigma_e, Es[::-1])
p = Seperate_Groups(Em, phi_vals, Es[::-1])[1]

E_log = e[::-1]
s     = s[::-1]
p     = p[::-1]

for g in range(len(Es)-1):
    E_g = E_log[g]
    for gp in range(len(Es)-1):
        E_gp = E_log[gp]
        vals = []
        for i in range(len(E_gp)):
            vals.append(np.trapz(f(E_g, alpha, E_gp[i]), E_g))
        phi_gp = np.trapz(p[gp], E_gp)
        R_gp_g = np.trapz(s[gp][i]*p[gp][i]*np.array(vals), E_gp)
        R[g, gp] = R_gp_g
        S[g, gp] = R[g, gp] / phi_gp 

plt.matshow(S)

a = S[0]
for i in range(1, len(S)):
    a += S[i]