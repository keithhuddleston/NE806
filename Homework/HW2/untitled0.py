# Imports ====================================================================
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# ============================================================================
A, mul, muc = sp.symbols('A \mu_L \mu_c')

expr = (A*muc+1) / sp.sqrt(A**2+2*A*muc+1)

a = sp.simplify(sp.diff(expr, muc)**-1)

b = sp.simplify(a.subs(A, 1))
