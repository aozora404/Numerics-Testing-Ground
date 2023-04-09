# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 14:52:17 2023

@author: Aozora

SI Units
"""

import matplotlib as mlib
import matplotlib.pyplot as plt
import numpy as np

"""
class Cell:
    
    position = np.array([0,0])
    
    electric_field = np.array([0, 0])
    magnetic_field = np.array(0)
    charge = 0
    current = np.array([0, 0])
    conductivity = 0

    def __init__(self, Position):
        self.position = Position

    def set_electric_field(self, Electric_field):
        self.electric_field = Electric_field

    def set_magnetic_field(self, Magnetic_field):
        self.magnetic_field = Magnetic_field

    def set_charge(self, Charge):
        self.charge = Charge

    def set_current(self, Current):
        self.current = Current
    
    def set_conductivity(self, Conductivity):
        self.conductivity = Conductivity
"""

c_squared = 299792458000 ** 2
epsilon_0_inverse = 1/(8.85418782 * 10**(-12))

grid_size = 1
cell_size = 0.1

cell_area = cell_size ** 2

cell_count = int(grid_size/cell_size)

# Initialize grid
electric_field = np.zeros(f"({cell_count}, {cell_count})")


        
for cells in grid:
    if np.sqrt(np.dot(cells.position, cells.position)) < 0.2:
        cells.conductivity = 20000
    else:
        cells.conductivity = 0

        