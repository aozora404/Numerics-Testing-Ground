# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 12:37:12 2023

@author: 12019044

Solver Functions
"""
import numpy as np


def gradient_flux(scalar_function, x, y, z):
    return 0.5 * np.array([scalar_function[z, y, x + 1] - scalar_function[z, y, x - 1], 
                           scalar_function[z, y + 1, x] - scalar_function[z, y - 1, x], 
                           scalar_function[z + 1, y, x] - scalar_function[z - 1, y, x]])
    
def divergence_flux(vector_function, x, y, z):
    return 0.5 * (vector_function[z, y, x + 1, 0] - vector_function[z, y, x - 1, 0] + 
                  vector_function[z, y + 1, x, 1] - vector_function[z, y - 1, x, 1] +
                  vector_function[z + 1, y, x, 2] - vector_function[z - 1, y, x, 2])
    
def curl_flux(vector_function, x, y, z):
    return 0.5 * np.array([vector_function[z, y + 1, x, 2] - vector_function[z, y - 1, x, 2]
                         - vector_function[z + 1, y, x, 1] + vector_function[z - 1, y, x, 1],
                           vector_function[z + 1, y, x, 0] - vector_function[z - 1, y, x, 0]
                         - vector_function[z, y, x + 1, 2] + vector_function[z, y, x - 1, 2],
                           vector_function[z, y, x + 1, 1] - vector_function[z, y, x - 1, 1]
                         - vector_function[z, y + 1, x, 0] + vector_function[z, y - 1, x, 0]])

def dot_product(vector_a, vector_b):
    return (vector_a[1] * vector_b[1] + vector_a[2] * vector_b[2] + vector_a[3] * vector_b[3])

def cross_product(vector_a, vector_b):
    return np.array(vector_a[2] * vector_b[3] - vector_a[3] * vector_b[2],
                    vector_a[3] * vector_b[1] - vector_a[1] * vector_b[3],
                    vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1])


def gradient_flux_2d(scalar_function, x, y):
    return 0.5 * np.array([scalar_function[y, x + 1] - scalar_function[y, x - 1], 
                           scalar_function[y + 1, x] - scalar_function[y - 1, x]])
    
def divergence_flux_2d(vector_function, x, y):
    return 0.5 * (vector_function[y, x + 1, 0] - vector_function[y, x - 1, 0] + 
                  vector_function[y + 1, x, 1] - vector_function[y - 1, x, 1])
    
def curl_flux_2d(vector_function, x, y):
    return 0.5 * (vector_function[y, x + 1, 1] - vector_function[y, x - 1, 1]
                - vector_function[y + 1, x, 0] + vector_function[y - 1, x, 0])