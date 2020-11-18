""" Homework for NE806, Neutronics
    
    Author: Keith Huddleston
     Email: kdhuddle@ksu.edu
"""

# ============================================================================
# Import statements
# ============================================================================
import numpy as np
import matplotlib.pyplot as plt


def fun(xi, Sig):
    return(-np.log(xi)/Sig)


Sig = 10
r = 1/Sig
n = 1000000
step = 5000
dist = [fun(np.random.rand(), Sig) for i in range(int(n/step))]
plt.hist(dist, 20), plt.show()


escape = [1 if i>r else 0 for i in dist]
print(np.mean(escape))

m = []
for i in range(step):
    dist = [fun(np.random.rand(), Sig) for i in range(int(n/step))]
    escape = [1 if i>r else 0 for i in dist]
    m.append(np.mean(escape))
    
plt.hist(m, 50), plt.axvline(np.mean(m), ls='--', c='k')