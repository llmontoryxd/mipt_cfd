# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 12:30:49 2020

@author: alexc
"""
import numpy as np
from numba import jit, njit
from matplotlib import pyplot as plt

# @jit
def solver(model, nx, nt, u0_fun, xl = 0., xr = 1., T = 1.):
    
    # B.c. are periodic
    
    # allocate arrays
    # solution (integral averages)
    u = np.zeros((nx, model.ne))
    
    # fluxes at all faces
    flux = np.zeros((nx + 1, model.ne))
    # reconstructed values
    # "0" - left value, "1" - right value
    ulr = np.zeros((nx+1, model.ne, 2))
    
    # right-hand side
    rhs = np.zeros((nx, model.ne))
    dx = (xr - xl) / nx
    dt = T / nt
    # cell centers
    xs = np.linspace(xl + dx/2, xr - dx/2, nx)
    
    # set initial condition
    u = u0_fun(xs)
    for it in range(nt):
        # Reconstruction
        for ix in range(nx):
            ulr[ix, :, 1]     = u[ix, :]
            ulr[ix + 1, :, 0] = u[ix, :]
        # Apply boundary conditions
        # periodic b.c. for simplicity
        ulr[0, :, 0] =  u[-1, :]
        ulr[-1, :, 1] = u[0, :]
        
        # Compute fluxes
        for jf in range(nx + 1):
            flux[jf, :] = model.flux_rs(ulr[jf, :, 0], ulr[jf, :, 1])
        
        # Compute rhs
        for ix in range(nx):
            rhs[ix, :] = -(flux[ix+1, :] - flux[ix, :]) / dx
        
        # Update u
        u = u + dt * rhs
    return u, xs, dx, dt
    
    

