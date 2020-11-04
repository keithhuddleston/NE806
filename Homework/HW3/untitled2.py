from scipy.integrate import quad
import numpy as np
sigH = 20
phi = lambda E: 1/E
alpha = 0.0
Es = np.logspace(-5, 7, 6)[::-1]
f = lambda E, alpha, Ep: 1/(Ep*(1-alpha)) if (alpha*Ep <= E <= Ep) else 0.0
R = np.zeros((len(Es)-1,len(Es)-1))
S = R*0
for g in range(len(Es)-1):
    for gp in range(len(Es)-1):
        E_gp = np.linspace(Es[gp+1], Es[gp], 1000)
        vals = []
        for i in range(len(E_gp)):
            vals.append(quad(lambda E: f(E, alpha, E_gp[i]), Es[g+1], Es[g])[0])
        #phi_gp = np.trapz(1/E_gp, E_gp) 
        # phi_gp = int(1/E,e0,e1) = ln(E0)-ln(E1)
        phi_gp = np.log(Es[gp]/Es[gp+1])
        R_gp_g = np.trapz(np.array(vals), E_gp)
        R[g, gp] = R_gp_g
        #print(phi_gp)
        S[g, gp] = R[g, gp] / phi_gp 
import matplotlib.pyplot as plt
R

