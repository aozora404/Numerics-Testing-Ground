# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 12:34:04 2023

Wave Equation
"""

import matplotlib as mlib
import matplotlib.pyplot as plt
import numpy as np
import numericslib as nlib
import math


# Pre-processing
c = 2
c_squared = c**2

grid_size = 10
cell_size = 0.04

cell_area = cell_size ** 2
cell_size_inv = 1/cell_size

cell_count = int(grid_size/cell_size)

displacement = np.zeros((cell_count + 2, cell_count + 2))
displacement_delta = np.zeros((cell_count + 2, cell_count + 2))
current = np.zeros((cell_count + 2, cell_count + 2, 2))

displacement_next = np.zeros((cell_count + 2, cell_count + 2))
displacement_delta_next = np.zeros((cell_count + 2, cell_count + 2))
current_next = np.zeros((cell_count + 2, cell_count + 2, 2))

# Set initial values
position = np.zeros((cell_count + 2, cell_count + 2, 2))

for y in range(1, cell_count):
    for x in range(1, cell_count):
        position[y,x,0] = -grid_size/2 + cell_size * x 
        position[y,x,1] = -grid_size/2 + cell_size * y

for y in range(1, cell_count):
    for x in range(1, cell_count):
        if np.sqrt(np.dot(position[y,x], position[y,x])) <= 0.5:
            displacement[y][x] = 0.05 * math.cos(math.pi * np.sqrt(np.dot(position[y,x], position[y,x])))

plt.figure(figsize=(10,10))
plt.pcolormesh(displacement, vmin=-0.1, vmax=0.1, cmap=plt.colormaps['seismic'])
plt.show()


# Solver
time_end = 10

t = 0
dt = 0.001

while t < time_end:
    for y in range(1, cell_count):
        for x in range(1, cell_count):
            displacement_next[y, x] = displacement[y, x] + dt * displacement_delta[y, x]
            displacement_delta_next[y, x] = displacement_delta[y, x] + dt * cell_size_inv * c_squared * nlib.divergence_flux_2d(current, x, y)
            current_next[y, x] = cell_size_inv * nlib.gradient_flux_2d(displacement, x, y)
    
    displacement = displacement_next
    displacement_delta = displacement_delta_next
    current = current_next
    
    t += dt
    
    print(t)
    
    plt.figure(figsize=(10,10))
    plt.pcolormesh(displacement, vmin=-0.1, vmax=0.1, cmap=plt.colormaps['seismic'])
    plt.show()

# Post-processing
plt.figure(figsize=(10,10))
plt.pcolormesh(displacement, vmin=-0.1, vmax=0.1, cmap=plt.colormaps['seismic'])
plt.show()

plt.figure(figsize=(10,10))
plt.quiver(position[::5,::5,0], position[::5,::5,1], current[::5,::5,0], current[::5,::5,1])
plt.show()