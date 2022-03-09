#!/usr/bin/env python
# coding: utf-8

# In[1]:


# -*- coding: utf-8 -*-
"""
Example: how to use mesh module
"""

import mesh as flowmod_mesh
from mesh import Mesh
import numpy as np

import time


path = 'mesh-cyl/'

mesh = Mesh()
t1 = time.time()
mesh.read_starcd(path)
t2 = time.time()

print('Time for mesh reading = {0:5.2e} s'.format(t2 - t1))


# Create random data and write it using 'write_tecplot' function 
# After that check that Tecplot opens this file and properly shows mesh
data = np.zeros((mesh.nc, 1))

for ic in range(mesh.nc) :
    data[ic, 0] = mesh.cell_volumes[ic] 
#    data[ic, 0] = np.random.rand() * ic

    
flowmod_mesh.write_tecplot(mesh, data, "volumes.tec", ["V"])
#flowmod_mesh.write_tecplot(mesh, data, "first.tec", ["random number * index of cell"])


# In[ ]:




