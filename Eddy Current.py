# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 06:48:18 2023

Eddy Current Solver
"""

import matplotlib as mlib
import matplotlib.pyplot as plt
import numpy as np
import numericslib as nlib

# Pre-processing
permeability = 1
mass = 1
mass_inverse = 1/mass

coil_current = 1
coil_radius = 1
coil_radius_squared = coil_radius**2
coil_origin = np.array([1,0])

initial_velocity = 1

grid_size = 2
cell_size = 0.05

cell_volume = cell_size ** 3
cell_size_inv = 1/cell_size

cell_count = int(grid_size/cell_size)

r = np.zeros((cell_count + 2, cell_count + 2, cell_count + 2))
for z in range(0, cell_count + 2):
    for y in range(0, cell_count + 2):
        for x in range(0, cell_count + 2):
            if(x!=0 & y!=0 & z!=0):
                r[z,y,x] = (x**2 + y**2 + z**2)**(-0.5)
            

def inverse_relative_distance(x_0, y_0, z_0, x, y, z):
    return r[abs(x_0-x), abs(y_0-y), abs(z_0-z)]

def A_field_embedded(position):
    s_vector = position[0:1] - coil_origin
    s = np.sqrt(s_vector.dot(s_vector))
    if s > coil_radius_squared:
        return 0.5 * permeability * coil_current * np.array([s_vector[1]*(coil_radius/s)**2, -s_vector[0]*(coil_radius/s)**2, 0])
    else:
        return 0.5 * permeability * coil_current * np.array([s_vector[1], -s_vector[0], 0])

def A_field_induced(x_0,y_0,z_0, current):
    integral = np.zeros(3)
    for z in range(1, cell_count):
        for y in range(1, cell_count):
            for x in range(1, cell_count):
                if(x!=x_0 & y!=y_0 & z!=z_0):
                    integral += current[z,y,x] * inverse_relative_distance(x_0, y_0, z_0, x, y, z)
    
    return permeability/(4 * np.pi) * integral * cell_volume

current = np.zeros((cell_count + 2, cell_count + 2, cell_count + 2, 3))
conductivity = np.zeros((cell_count + 2, cell_count + 2, cell_count + 2))

magnetic_field = np.zeros((cell_count + 2, cell_count + 2, cell_count + 2, 3))

magnetic_potential = np.zeros((cell_count + 2, cell_count + 2, cell_count + 2, 3))
magnetic_potential_next = np.zeros((cell_count + 2, cell_count + 2, cell_count + 2, 3))
magnetic_potential_previous = np.zeros((cell_count + 2, cell_count + 2, cell_count + 2, 3))

force = np.zeros(3)
velocity = np.zeros(3)
velocity_next = np.zeros(3)
position = np.zeros(3)
position_next = np.zeros(3)

# Set initial values
velocity = np.array([initial_velocity, 0, 0])
relative_position = np.zeros((cell_count + 2, cell_count + 2, cell_count + 2, 3))

for z in range(1, cell_count):
    for y in range(1, cell_count):
        for x in range(1, cell_count):
            relative_position[z,y,x,0] = -grid_size/2 + cell_size * x 
            relative_position[z,y,x,1] = -grid_size/2 + cell_size * y
            relative_position[z,y,x,2] = -grid_size/2 + cell_size * z

for z in range(1, cell_count):
    for y in range(1, cell_count):
        for x in range(1, cell_count):
            if np.sqrt(np.dot(relative_position[z,y,x], relative_position[z,y,x])) <= 0.5:
                conductivity[z,y,x] = 1

for z in range(1, cell_count):
    for y in range(1, cell_count):
        for x in range(1, cell_count):
            magnetic_potential[z,y,x] = A_field_embedded(position + relative_position[z,y,x])
            magnetic_potential_previous[z,y,x] = A_field_embedded(position + relative_position[z,y,x])

plt.figure(figsize=(10,10))
plt.quiver(relative_position[::1,int(cell_count/2),::1,0], relative_position[::1, int(cell_count/2),::1,1], magnetic_potential[::1,int(cell_count/2),::1,0], magnetic_potential[::1,int(cell_count/2),::1,1])
plt.show()


# Solver
time_end = 1

t = 0
dt = 0.001

while t < time_end:
    for z in range(1, cell_count):
        for y in range(1, cell_count):
            for x in range(1, cell_count):
                current[z,y,x] = - conductivity[z,y,x] * (magnetic_potential[z,y,x] - magnetic_potential_previous[z,y,x]) / dt
                magnetic_field[z,y,x] = cell_size_inv * nlib.curl_flux(magnetic_potential, x, y, z)
    
    force = np.zeros(3)
    for z in range(1, cell_count):
        for y in range(1, cell_count):
            for x in range(1, cell_count):
                force += np.cross(current[z,y,x], magnetic_field[z,y,x])
    
    force *= cell_volume
    velocity_next = velocity + dt * force * mass_inverse
    position_next = position + dt * velocity 
    
    for z in range(1, cell_count):
        for y in range(1, cell_count):
            for x in range(1, cell_count):
                magnetic_potential_next[z,y,x] = A_field_embedded(position + relative_position[z,y,x]) + A_field_induced(x,y,z, current)
    
    magnetic_potential_previous = magnetic_potential
    magnetic_potential = magnetic_potential_next
    
    
    
    t += dt
    
    print(t)
    
    plt.figure(figsize=(10,10))
    plt.quiver(relative_position[::5,int(cell_count/2),::5,0], relative_position[::5,int(cell_count/2),::5,1], current[::5,int(cell_count/2),::5,0], current[::5,int(cell_count/2),::5,1])
    plt.show()
    

# Post-processing
plt.figure(figsize=(10,10))
plt.quiver(relative_position[::5,int(cell_count/2),::5,0], relative_position[::5,int(cell_count/2),::5,1], magnetic_potential[::5,int(cell_count/2),::5,0], magnetic_potential[::5,int(cell_count/2),::5,1])
plt.show()

plt.figure(figsize=(10,10))
plt.quiver(relative_position[::5,int(cell_count/2),::5,0], relative_position[::5,int(cell_count/2),::5,1], current[::5,int(cell_count/2),::5,0], current[::5,int(cell_count/2),::5,1])
plt.show()