# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 12:37:12 2023

@author: 12019044

Solver Functions
"""
import numpy as np


def gradient_flux(scalar_function, x, y, z, cell_area):
    return cell_area * 0.5 * np.array([scalar_function[z, y, x + 1] - scalar_function[z, y, x - 1], 
                                       scalar_function[z, y + 1, x] - scalar_function[z, y - 1, x], 
                                       scalar_function[z + 1, y, x] - scalar_function[z - 1, y, x]])
    
def divergence_flux(vector_function, x, y, z, cell_area):
    return cell_area * 0.5 * (vector_function[z, y, x + 1, 0] - vector_function[z, y, x - 1, 0] + 
                              vector_function[z, y + 1, x, 1] - vector_function[z, y - 1, x, 1] +
                              vector_function[z + 1, y, x, 2] - vector_function[z - 1, y, x, 2])
    
def curl_flux(vector_function, x, y, z, cell_area):
    return cell_area * 0.5 * np.array([vector_function[z, y + 1, x, 2] - vector_function[z, y - 1, x, 2]
                                     - vector_function[z + 1, y, x, 1] + vector_function[z - 1, y, x, 1],
                                       vector_function[z + 1, y, x, 0] - vector_function[z - 1, y, x, 0]
                                     - vector_function[z, y, x + 1, 2] + vector_function[z, y, x - 1, 2],
                                       vector_function[z, y, x + 1, 1] - vector_function[z, y, x - 1, 1]
                                     - vector_function[z, y + 1, x, 0] + vector_function[z, y - 1, x, 0]])


def gradient_flux_2d(scalar_function, x, y, cell_size):
    return cell_size * 0.5 * np.array([scalar_function[y, x + 1] - scalar_function[y, x - 1], 
                                       scalar_function[y + 1, x] - scalar_function[y - 1, x]])
    
def divergence_flux_2d(vector_function, x, y, cell_size):
    return cell_size * 0.5 * (vector_function[y, x + 1, 0] - vector_function[y, x - 1, 0] + 
                              vector_function[y + 1, x, 1] - vector_function[y - 1, x, 1])
    
def curl_flux_2d(vector_function, x, y, cell_size):
    return cell_size * 0.5 * (vector_function[y, x + 1, 1] - vector_function[y, x - 1, 1]
                            - vector_function[y + 1, x, 0] + vector_function[y - 1, x, 0])