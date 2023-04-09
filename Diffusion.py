# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 06:48:18 2023

Diffusion Equation Solver
"""

import matplotlib as mlib
import matplotlib.pyplot as plt
import numpy as np


# Pre-processing
alpha = 100000000

grid_size = 10
cell_size = 0.05

cell_area = cell_size ** 2

cell_count = int(grid_size/cell_size)

charge = np.zeros((cell_count + 2, cell_count + 2))
current = np.zeros((cell_count + 2, cell_count + 2, 2))
charge_next = np.zeros((cell_count + 2, cell_count + 2))
current_next = np.zeros((cell_count + 2, cell_count + 2, 2))

# Set initial values
position = np.zeros((cell_count + 2, cell_count + 2, 2))

for y in range(1, cell_count):
    for x in range(1, cell_count):
        position[y,x,0] = -grid_size/2 + cell_size * x 
        position[y,x,1] = -grid_size/2 + cell_size * y

for y in range(1, cell_count):
    for x in range(1, cell_count):
        if np.sqrt(np.dot(position[y,x], position[y,x])) <= 2:
            charge[y][x] = 0.05

plt.figure(figsize=(10,10))
plt.pcolormesh(charge)
plt.show()


# Solver
time_end = 30

t = 0
dt = 0.1

while t < time_end:
    for y in range(1, cell_count):
        for x in range(1, cell_count):
            charge_next[y, x] = charge[y, x] + dt * (-cell_area * cell_size * 0.5 * (current[y, x + 1, 0] + current[y + 1, x, 1] - current[y, x - 1, 0] - current[y - 1, x, 1]))
            
            current_next[y, x] = -alpha * cell_area * cell_size * 0.5 * np.array([charge[y, x + 1] - charge[y, x - 1], charge[y+1, x] - charge[y-1, x]])
    
    charge = charge_next
    current = current_next
    
    t += dt
    
    print(t)
    

# Post-processing
plt.figure(figsize=(10,10))
plt.pcolormesh(charge)
plt.show()

plt.figure(figsize=(10,10))
plt.quiver(position[::5,::5,0], position[::5,::5,1], current[::5,::5,0], current[::5,::5,1])
plt.show()