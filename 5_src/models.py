# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 11:56:40 2020

@author: alexc
"""

import numpy as np


# Linear advection equation: u_t + a u_x = 0
class Advection():
    def __init__(self, a = 1.):
        if (a == 0.):
            raise NameError('a equals 0!')
        self.a = a
        self.ne = 1
    
    # Flux function
    def flux(self, u):
        return self.a * u
    
    # Riemann problem solver
    def flux_rs(self, ul, ur):
        """

        Parameters
        ----------
        ul : np.array
            vectors of unknows to the left of a cell face
        ur : np.array
            vectors of unknows to the right of a cell face

        Returns
        -------
        np.array
            fluxes at faces computed via solution of the Riemann problem

        """
        if (self.a > 0):
            return self.flux(ul)
        else:
            return self.flux(ur)
        
# Inviscid Burgers equation: u_t + (u^2 / 2)_x = 0
# TODO: add your code

