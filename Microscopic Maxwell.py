# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 12:19:41 2023

@author: 12019044
"""

import matplotlib as mlib
import matplotlib.pyplot as plt
import numpy as np


# Pre-processing
alpha = 100000000

grid_size = 10
cell_size = 0.05

cell_volume = cell_size ** 3

cell_count = int(grid_size/cell_size)

e_field = np.zeros((cell_count + 2, cell_count + 2, 3))
e_field_delta = np.zeros((cell_count + 2, cell_count + 2, 3))
e_field_divergence = np.zeros((cell_count + 2, cell_count + 2))
e_field_curl = np.zeros((cell_count + 2, cell_count + 2, 3))

b_field = np.zeros((cell_count + 2, cell_count + 2, 3))
b_field_delta = np.zeros((cell_count + 2, cell_count + 2, 3))
b_field_divergence = np.zeros((cell_count + 2, cell_count + 2))
b_field_curl = np.zeros((cell_count + 2, cell_count + 2, 3))

charge = np.zeros((cell_count + 2, cell_count + 2))
current = np.zeros((cell_count + 2, cell_count + 2, 3))


e_field_next = np.zeros((cell_count + 2, cell_count + 2, 3))
e_field_delta_next = np.zeros((cell_count + 2, cell_count + 2, 3))
e_field_divergence_next = np.zeros((cell_count + 2, cell_count + 2))
e_field_curl_next = np.zeros((cell_count + 2, cell_count + 2, 3))

b_field_next = np.zeros((cell_count + 2, cell_count + 2, 3))
b_field_delta_next = np.zeros((cell_count + 2, cell_count + 2, 3))
b_field_divergence_next = np.zeros((cell_count + 2, cell_count + 2))
b_field_curl_next = np.zeros((cell_count + 2, cell_count + 2, 3))

charge_next = np.zeros((cell_count + 2, cell_count + 2))
current_next = np.zeros((cell_count + 2, cell_count + 2, 3))


# Set initial values


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
plt.pcolormesh(charge, vmin=0, vmax=0.05)
plt.show()

plt.figure(figsize=(10,10))
plt.quiver(position[::5,::5,0], position[::5,::5,1], current[::5,::5,0], current[::5,::5,1])
plt.show()