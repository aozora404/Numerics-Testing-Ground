# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 17:38:14 2023

@author: Aozora
"""

import numpy as np
import matplotlib.pyplot as plt

# Generate some test data
data = np.arange(100).reshape((10,10))

plt.title('Actual Function')
heatmap = plt.pcolormesh(data)
plt.show()