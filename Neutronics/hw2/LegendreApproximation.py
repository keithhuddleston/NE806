"""
File written for NE 806, Neutronics HW2
"""


#  Imports  ==================================================================
import numpy as np
import sympy as sp
sp.init_printing(use_unicode=True)
import matplotlib.pyplot as plt


#  Define functions for legendre approximation
def legendrePol(x, n):
    return (1/(2**n*sp.factorial(n))) * sp.diff((x**2-1)**n, x, n)

def legendreA(x, n, function, LP):
    return ((2*n + 1)/2) * sp.integrate(function*LP, (x, -1, 1))

def approx(degree):
    approximation = 0
    for i in range(degree):
        LP = legendrePol(muc, i)
        approximation += LP*legendreA(muc, i, b, LP)
    return approximation


#  Define P(\mu_L)  ==========================================================
A, mul, muc = sp.symbols('A \mu_L \mu_c')

expr = (A*muc+1) / sp.sqrt(A**2+2*A*muc+1)

a = sp.simplify(sp.diff(expr, muc)**-1)

b = sp.simplify(a.subs(A, 1))/4/sp.pi

# First four legendre approximations
A0 = approx(0)
A1 = approx(1)
A2 = approx(2)
A3 = approx(3)

# Convert sympy functions to functions which may be evaluated
PA0 = np.vectorize(sp.lambdify(muc, A0, "numpy"))
PA1 = np.vectorize(sp.lambdify(muc, A1, "numpy"))
PA2 = np.vectorize(sp.lambdify(muc, A2, "numpy"))
PA3 = np.vectorize(sp.lambdify(muc, A3, "numpy"))
base = sp.lambdify(muc, b, "numpy") # Equation that we are approximated

# Plotting ===================================================================
plt.rcParams.update({'font.size': 18})
z = np.linspace(-1, 1, 200)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))

ax.plot(z, base(z), z, PA0(z), z, PA1(z), z, PA2(z), z, PA3(z))

plt.tight_layout(), plt.grid()
plt.xlabel('$\mu_c$'), plt.ylabel('$P(\mu_L)$')
plt.legend(['Exact Function', '$P_0$', '$P_1$', '$P_2$', '$P_3$'])

