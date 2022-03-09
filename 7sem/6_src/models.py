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
        

# 1D Euler (Gas Dynamics) equations
class Euler():
    def __init__(self, rs_type = 'hll', wave_speed_estimate_type = 'Davis'):
        self.ne = 3
        self.gamma = 1.4
        # Which Riemann solver is used
        self.rs_type = rs_type
        self.wave_speed_estimate_type = wave_speed_estimate_type
    
    # compute primitive variables V[:]
    # from vector of conservative variables U[:]   
    def cons2prim(self, U):
        V = np.zeros_like(U)
        V[0] = U[0] # rho
        V[1] = U[1] / U[0] # u = rho * u / u
        # p = (E - 0.5 rho u^2) (gamma - 1)
        V[2] = (U[2] - 0.5 * V[0] * V[1]**2) * (self.gamma - 1)
        return V
        
    # compute  conservative variables U[:]
    # from vector of primitive variables V[:]   
    def prim2cons(self, V):
        U = np.zeros_like(V)
        U[0] = V[0] # rho
        U[1] = V[0] * V[1] # rho * u
        # E = p/(gamma - 1) + 0.5 rho u^2
        U[2] = V[2] / (self.gamma - 1) + 0.5 * V[0] * V[1]**2
        return U
    
    # Flux function
    def flux(self, U):
        # U - vector of conservative variables
        pass
    
    def flux_hll(self, UL, UR):
        pass
        
    # Riemann problem solver
    def flux_rs(self, UL, UR):
        if(self.rs_type =='hll'):
            return self.flux_hll(UL, UR)
    
    def speed_of_sound(self, U):
        V = self.cons2prim(U)
        # (gamma * p / rho)^1/2
        return np.sqrt(self.gamma * V[2] / V[0])
    
    def davis_estimate(self, UL, UR):
        VL = self.cons2prim(UL)
        VR = self.cons2prim(UR)
        
        aL = self.speed_of_sound(UL)
        aR = self.speed_of_sound(UR)
        SL = min(VL[1] - aL, VR[1] - aR)
        SR = max(VL[1] + aL, VR[1] + aR)
        return SL, SR
    
    def einfeldt_estimate(self, UL, UR):    
        pass
        
    def wave_speed_estimate(self, UL, UR):
        if (self.wave_speed_estimate_type == 'Davis'):
            return self.davis_estimate(UL, UR)
        elif (self.wave_speed_estimate_type == 'Einfeldt'):
            return self.einfeldt_estimate(UL, UR)
        
        
