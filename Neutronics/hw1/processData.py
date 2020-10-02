# Module Imports =============================================================
import numpy as np
import matplotlib.pyplot as plt
import re

file = open('u235.txt', 'r')

lines = file.readlines()
l = len(lines)  # Number of Rows
c = 6  # Number of Columns
data = np.zeros((l, c))

p1 = r'(\+*\-*\d\.\d*\+*\-*\d+)'
for i in range(len(lines)):
    line = lines[i]
    found = re.findall(p1, line)
    for j in range(c):
        value = found[j]
        for k in range(len(value))[::-1]:
            if value[k] == '-' or value[k] == '+':
                value = value[:k]+'e'+value[k:]
                break
        data[i][j] = value

np.savetxt("u235processed.txt", data)
    