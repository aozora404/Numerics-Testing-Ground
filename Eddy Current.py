# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 06:48:18 2023

Eddy Current Solver
"""

import matplotlib as mlib
import matplotlib.pyplot as plt
import numpy as np
import numericslib as nlib

# Settings
# Simulation settings
grid_size = 22 # mm
cell_size = 0.088 # mm
time_step = 0.000008 # s
time_end = 1 # s

# Flying mass settings
initial_velocity = 0 # mm s-1
mass = 0.01 # kg
object_radius = 10 # mm 
material_conductivity = 3.5 * 10**1 # mm-2 kg-1 s3 A2
material_magnetic_permeability = 1.256665 * 10**(-3) # mm kg s-2 A-2

# Magnetic coil settings
coil_radius = 200 # mm
coil_magnetic_field_strength = -0.1 # kg s-2 A-1
coil_distance = 208 # mm

# Miscellaneous
force_external_magnitude = 5000 # mm kg s-2


# Pre-processing
# Constants and Coefficients
mu_0 = 1.25663706212 * 10**(-3) # kg mm s-2 A-2
mu_0_inverse = 1/mu_0

permeability = material_magnetic_permeability
permeability_inverse = 1/permeability

mass_inverse = 1/mass

coil_radius_squared = coil_radius**2
coil_origin = np.array([coil_distance,0])

alpha = 0.5 * coil_magnetic_field_strength


cell_volume = cell_size ** 2
cell_size_inv = 1/cell_size

cell_count = int(grid_size/cell_size)


# Initialize field grid
current = np.zeros((cell_count + 2, cell_count + 2, 2))
conductivity = np.zeros((cell_count + 2, cell_count + 2))

magnetic_field = np.zeros((cell_count + 2, cell_count + 2))


magnetic_potential = np.zeros((cell_count + 2, cell_count + 2, 2))
magnetic_potential_next = np.zeros((cell_count + 2, cell_count + 2, 2))

magnetic_potential_current_x = np.zeros((cell_count + 2, cell_count + 2, 2))
magnetic_potential_current_x_next = np.zeros((cell_count + 2, cell_count + 2, 2))

magnetic_potential_current_y = np.zeros((cell_count + 2, cell_count + 2, 2))
magnetic_potential_current_y_next = np.zeros((cell_count + 2, cell_count + 2, 2))

magnetic_potential_delta = np.zeros((cell_count + 2, cell_count + 2, 2))
magnetic_potential_delta_next = np.zeros((cell_count + 2, cell_count + 2, 2))

magnetic_potential_embedded = np.zeros((cell_count + 2, cell_count + 2, 2))
magnetic_potential_embedded_previous = np.zeros((cell_count + 2, cell_count + 2, 2))

magnetic_potential_embedded_delta = np.zeros((cell_count + 2, cell_count + 2, 2))
magnetic_potential_embedded_lagrangian = np.zeros((cell_count + 2, cell_count + 2, 2))


force = np.zeros(2)
force_external = np.zeros(2)
velocity = np.zeros(2)
velocity_next = np.zeros(2)
position = np.zeros(2)
position_next = np.zeros(2)

graph_time = np.array([])
graph_position = np.array([])
graph_velocity = np.array([])
graph_force = np.array([])

# Set initial values
force_external = np.array([force_external_magnitude, 0])
velocity = np.array([initial_velocity, 0])
relative_position = np.zeros((cell_count + 2, cell_count + 2, 2))

for y in range(1, cell_count):
    for x in range(1, cell_count):
        relative_position[y,x,0] = -grid_size/2 + cell_size * x
        relative_position[y,x,1] = -grid_size/2 + cell_size * y

for y in range(1, cell_count):
    for x in range(1, cell_count):
        if np.dot(relative_position[y,x], relative_position[y,x]) <= object_radius**2:
            conductivity[y,x] = material_conductivity

for y in range(1, cell_count):
    for x in range(1, cell_count):
        s = relative_position[y,x] + position - coil_origin
        s_mag = s.dot(s)
        if s_mag > coil_radius_squared:
            magnetic_potential_embedded[y,x] = alpha * (coil_radius_squared/s_mag) * np.array([-s[1], s[0]])
        else:
            magnetic_potential_embedded[y,x] = alpha * np.array([-s[1], s[0]])

magnetic_potential_embedded_previous = magnetic_potential_embedded


plt.figure(figsize=(10,10))
plt.quiver(relative_position[::4,::4,0], relative_position[::4,::4,1], magnetic_potential[::4,::4,0]+magnetic_potential_embedded[::4,::4,0], magnetic_potential[::4,::4,1]+magnetic_potential_embedded[::4,::4,1])
plt.show()


# Solver
t = 0
dt = time_step

while t < time_end:
    
    # Calculate embedded magnetic potential
    for y in range(1, cell_count):
        for x in range(1, cell_count):
            s = relative_position[y,x] + position - coil_origin
            s_mag = s.dot(s)
            if s_mag > coil_radius_squared:
                magnetic_potential_embedded[y,x] = alpha * (coil_radius_squared/s_mag) * np.array([-s[1], s[0]])
                magnetic_potential_embedded_delta[y,x] = (1/dt)*(magnetic_potential_embedded[y,x] - magnetic_potential_embedded_previous[y,x])
                magnetic_potential_embedded_lagrangian[y,x] = magnetic_potential_embedded[y,x] / s_mag
            else:
                magnetic_potential_embedded[y,x] = alpha * np.array([-s[1], s[0]])
                magnetic_potential_embedded_delta[y,x] = alpha * np.array([-velocity[1], velocity[0]])
                magnetic_potential_embedded_lagrangian[y,x] = magnetic_potential_embedded[y,x] / s_mag
            
    # Calculate induced magnetic potential
    for y in range(1, cell_count):
        for x in range(1, cell_count):
            if conductivity[y,x] != 0:
                magnetic_potential_delta_next[y,x] = (permeability_inverse/conductivity[y,x]) * (cell_size_inv * np.array([nlib.divergence_flux_2d(magnetic_potential_current_x, x, y), nlib.divergence_flux_2d(magnetic_potential_current_y, x, y)]) + magnetic_potential_embedded_lagrangian[y,x]) - magnetic_potential_embedded_delta[y,x] 
                magnetic_potential_current_x_next[y,x] = cell_size_inv * nlib.gradient_flux_2d(magnetic_potential[:,:,0], x, y)
                magnetic_potential_current_y_next[y,x] = cell_size_inv * nlib.gradient_flux_2d(magnetic_potential[:,:,1], x, y)
                magnetic_potential_next[y,x] = magnetic_potential[y,x] + dt * magnetic_potential_delta[y,x]

    # Current and Field
    for y in range(1, cell_count):
        for x in range(1, cell_count):
            if conductivity[y,x] != 0:
                current[y,x] = - conductivity[y,x] * (magnetic_potential_delta[y,x] + magnetic_potential_embedded_delta[y,x])
                magnetic_field[y,x] = cell_size_inv * nlib.curl_flux_2d(magnetic_potential + magnetic_potential_embedded , x, y)
    
    # Dynamics
    force = np.zeros(2)
    for y in range(1, cell_count):
        for x in range(1, cell_count):
            if conductivity[y,x] != 0:
                force += magnetic_field[y,x] * np.array([current[y,x,1], -current[y,x,0]])
    
    force *= cell_volume
    
    force += force_external
    
    velocity_next = velocity + dt * force * mass_inverse
    position_next = position + dt * velocity 
    
    # Time step
    magnetic_potential_previous = magnetic_potential
    magnetic_potential = magnetic_potential_next
    magnetic_potential_current_x = magnetic_potential_current_x_next
    magnetic_potential_current_y = magnetic_potential_current_y_next
    magnetic_potential_delta = magnetic_potential_delta_next
    magnetic_potential_embedded_previous = magnetic_potential_embedded
    position = position_next
    velocity = velocity_next

    t += dt
    
    graph_time = np.append(graph_time, [t])
    graph_position = np.append(graph_position, [position[0]])
    graph_velocity = np.append(graph_velocity, [velocity[0]])
    graph_force = np.append(graph_force, [force[0]])
    
    print(t)
    print("Position = ", position)
    print("Velocity = ", velocity)
    print("Force = ", force)
    print("")
    
    plt.figure(figsize=(10,10))
    plt.quiver(relative_position[::4,::4,0], relative_position[::4,::4,1], current[::4,::4,0], current[::4,::4,1])
    plt.show()
    """
    plt.figure(figsize=(10,10))
    plt.pcolormesh(magnetic_field, cmap=plt.colormaps['seismic'], vmin=-0.001, vmax=0.001)
    plt.show()
    """

# Post-processing
plt.figure(figsize=(10,10))
plt.scatter(graph_time, graph_position)
plt.show()

plt.figure(figsize=(10,10))
plt.scatter(graph_time, graph_velocity)
plt.show()

plt.figure(figsize=(10,10))
plt.scatter(graph_time, graph_force)
plt.show()

plt.figure(figsize=(10,10))
plt.quiver(relative_position[::4,::4,0], relative_position[::4,::4,1], magnetic_potential[::4,::4,0], magnetic_potential[::4,::4,1])
plt.show()

plt.figure(figsize=(10,10))
plt.quiver(relative_position[::4,::4,0], relative_position[::4,::4,1], current[::4,::4,0], current[::4,::4,1])
plt.show()