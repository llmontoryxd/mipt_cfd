# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 11:56:40 2020

@author: alexc
"""

import numpy as np

# ---------------------------------------------------------------------------- #
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

# ---------------------------------------------------------------------------- #
# Inviscid Burgers equation: u_t + (u^2 / 2)_x = 0
class Burgers() :
    def __init__(self) :
        self.ne = 1
        
    # Flux function
    def flux(self, u) :
        return u**2/2
    
    # Riemann problem solver
    def flux_rs(self, ul, ur) :
        urs = False
        
        if (ul[0] > ur[0]) :
            S = 0.5 * (ul + ur)
            if (S[0] > 0) :
                urs = ul
            else :
                urs = ur
        else :
            if (ul[0] >= 0) :
                urs = ul
            elif (ur[0] <= 0) :
                urs = ur
            else :
                urs = 0
                
        if (urs == False) :
            print("Some problem with Riemann problem solver!")
            
        return self.flux(urs)

# ---------------------------------------------------------------------------- #
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
    
    # flux function
    def flux(self, U):
        # U - vector of conservative variables
        # V - vector of primitive variables
        V = self.cons2prim(U)
        # f - vector of flux
        f = np.zeros_like(U)

        # compute flux 
        f[0] = U[1]
        f[1] = U[1] * V[1] + V[2]
        f[2] = V[1] * (U[2] + V[2])

        return f
        
    # flux function in HLL method
    def flux_hll(self, UL, UR):
        # SL and SR - min and max speed of wave
        SL, SR = self.wave_speed_estimate(UL, UR)
        
        if (SL > 0) :
            return self.flux(UL)
        elif (SR < 0) : 
            return self.flux(UR)
        else :
            # FL = flux(UL), FR = flux(UR)
            FL = self.flux(UL)
            FR = self.flux(UR)
            return (SR * FL - SL * FR  + SL * SR * (UR - UL)) / (SR - SL)

    # flux function in HLLC method
    def flux_hllc(self, UL, UR) : 
        # SL and SR - min and max speed of wave
        SL, SR = self.wave_speed_estimate(UL, UR)
        
        # SS - S* - intermediate wave
        SS = self.SS(UL, UR)

        if (SL >= 0) :
            return self.flux(UL)
        elif (SS >= 0 and SL <= 0) : 
            USL = self.US(UL, SL, SS)
        
            return self.FS(SL, USL, UL)
        elif (SR >= 0 and SS <= 0) : 
            USR = self.US(UR, SR, SS)

            return self.FS(SR, USR, UR)
        else :
            return self.flux(UR)

    # S* function
    def SS(self, UL, UR) :
        # VL and VR - vectors of primitive variables from UL and UR
        VL = self.cons2prim(UL)
        VR = self.cons2prim(UR)

        # SL and SR - min and max speed of wave
        SL, SR = self.wave_speed_estimate(UL, UR)

        # here pR = VR[2], pL = VL[2], rhoL*uL = UL[1], rhoR*uR = UR[1], uL = VL[1], uR = VR[1], rhoL = VL[0], rhoR = VR[0]
        return (VR[2] - VL[2] + UL[1] * (SL - VL[1]) - UR[1] * (SR - VR[1])) / (VL[0] * (SL - VL[1]) - VR[0] * (SR - VR[1]))

    # U* function        
    def US(self, U, S, SS) : 
        # V - vector of primitive variables
        V = self.cons2prim(U)

        # here rho = V[0], u = V[1]
        coef = V[0] * (S - V[1]) / (S - SS)

        # U in star region
        US = np.ones_like(U) * coef

        # here E = U[2], p = V[2]
        US[1] = US[1] * SS
        US[2] = US[2] * (U[2]/V[0] + (SS - V[1]) * (SS + V[2]/(V[0] * (S - V[1]))))

        return US

    # F* function
    def FS(self, S, US, U) : 
        return self.flux(U) + S*(US - U)

    # Riemann problem solver
    def flux_rs(self, UL, UR):
        if(self.rs_type == 'hll'):
            return self.flux_hll(UL, UR)
        if (self.rs_type == 'hllc'):
            return self.flux_hllc(UL, UR)
    
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
        # VL and VR - vectors of primitive variables from UL and UR 
        VL = self.cons2prim(UL)
        VR = self.cons2prim(UR)

        # u - the ROE-averaged value of u
        u = (np.sqrt(VL[0]) * VL[1] + np.sqrt(VR[0]) * VR[1]) / (np.sqrt(VL[0]) + np.sqrt(VR[0]))
        
        # Enthalpy on the left and right
        HL = (UL[2] + VL[2]) / VL[0]
        HR = (UR[2] + VR[2]) / VR[0]
        
        # H - the ROE-averaged value of H
        H = (np.sqrt(VL[0]) * HL + np.sqrt(VR[0]) * HR) / (np.sqrt(VL[0]) + np.sqrt(VR[0]))

        # a - the ROE-averaged value of a
        a = np.sqrt((self.gamma - 1) * (H - 0.5 * u**2))

        return u-a, u+a

        
    def wave_speed_estimate(self, UL, UR):
        if (self.wave_speed_estimate_type == 'Davis'):
            return self.davis_estimate(UL, UR)
        elif (self.wave_speed_estimate_type == 'Einfeldt'):
            return self.einfeldt_estimate(UL, UR)
        
        
