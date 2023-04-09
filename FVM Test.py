# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 15:29:56 2023

@author: Aozora
"""

"""
'''solves''' f' + g = 0
"""

from math import *

f = []

def g(x):
    return -1


dx = 0.0001
x = 2

n = int(x/dx)


for i in range(n):
    f.append(0)

f[0] = 0


for i in range(n-1):
    f[i+1] = f[i] - g(i)*dx
    
print(f[n-1])