# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 11:43:59 2020

@author: alexc
"""
import models
from solver import solver
import numpy as np
from matplotlib import pyplot as plt
import time

################################
# Test for advection equation
################################

# Create model for advection equation
adv_model = models.Advection()

# Define function for initial conditions
def u0_1(x):
    u0 = np.zeros((len(x), 1))
    u0[:, 0] = np.exp(-(x - 0.5)**2 / 0.01)
    return u0

# Plot initial condition
xs = np.linspace(0, 1, 50)
plt.plot(xs, u0_1(xs), label = 'exact')

# Compute numerical solution
t1 = time.time()
u, xs = solver(adv_model, nx = 50, nt=100, T = 1, u0_fun = u0_1)
t2 = time.time()
print ('Solution computed in {0:5.2e} seconds'.format(t2 - t1))
# Plot numerical solution
plt.plot(xs, u[:, 0], 'o', label = 'Numerical')
plt.grid(True)
plt.legend()
plt.show()

################################
# Test for Burgers equation
################################

# Create model for advection equation
burg_model = models.Burgers()

# Define function for initial conditions
def u0_2(x):
    u0 = np.zeros((len(x), 1))
    u0[:, 0] = 1 + 0.5 * np.sin(np.pi * x)
    return u0

# Compute numerical solution
t1 = time.time()
u, xs = solver(burg_model, nx = 50, nt=100, T = 1.1, u0_fun = u0_2,
               xl = -1., xr = + 1.)
t2 = time.time()
print ('Solution computed in {0:5.2e} seconds'.format(t2 - t1))
# Plot numerical solution
plt.plot(xs, u[:, 0], '-o', label = 'Numerical, T = 1.1')
plt.grid(True)
plt.legend()
# plt.savefig('burgers_t_1_1.png', dpi = 600, bbox_inches = 'tight')
plt.show()

###############################
# Tests for Euler equations
###############################

##################################
## Test 2
##################################

euler_model = models.Euler(wave_speed_estimate_type = 'Einfeldt')

# Define left and right states for the test Riemann problem
VL = np.array([1., 0.75, 1.0])
VR = np.array([0.125, 0.0, 0.1])

# Point of discontinuity 
x0 = 0.3
T_final = 0.2

UL = euler_model.prim2cons(VL)
UR = euler_model.prim2cons(VR)
# Define function for initial conditions
def u0_euler(x):
    U0 = np.zeros((len(x), 3))
    mask = x <= x0
    U0[mask, :] = UL
    mask = x > x0
    U0[mask, :] = UR

    return U0


xx = np.linspace(0, 1, 200)
U0 = u0_euler(xx)
# print("U0.shape = ", U0.shape)
# plt.plot(xx, U0[:, 0])

t1 = time.time()
U, xs = solver(euler_model, nx = 100, nt=800, T = T_final, u0_fun = u0_euler,
               xl = 0., xr = + 1., bc_type = 'transparent')
t2 = time.time()
print ('Solution computed in {0:5.2e} seconds'.format(t2 - t1))

plt.plot(xs, U[:, 0], 'o-')
plt.grid(True)
plt.show()


##################################
## Test 2
##################################
# euler_model = models.Euler(wave_speed_estimate_type = 'Davis')
# # Define left and right states for the test Riemann problem
# VL = np.array([1., -2, 0.4])
# VR = np.array([1, 2.0, 0.4])

# # Point of discontinuity 
# x0 = 0.5
# T_final = 0.15

# UL = euler_model.prim2cons(VL)
# UR = euler_model.prim2cons(VR)
# # Define function for initial conditions
# def u0_euler(x):
#     U0 = np.zeros((len(x), 3))
#     mask = x <= x0
#     U0[mask, :] = UL
#     mask = x > x0
#     U0[mask, :] = UR

#     return U0


# xx = np.linspace(0, 1, 200)
# U0 = u0_euler(xx)
# print("U0.shape = ", U0.shape)
# # plt.plot(xx, U0[:, 0])

# t1 = time.time()
# U, xs = solver(euler_model, nx = 100, nt=400, T = T_final, u0_fun = u0_euler,
#                 xl = 0., xr = + 1., bc_type = 'transparent')
# t2 = time.time()
# print ('Solution computed in {0:5.2e} seconds'.format(t2 - t1))

# plt.plot(xs, U[:, 0], 'o-')
# plt.grid(True)
# plt.show()


# ##################################
# ## Test 3
# ##################################
# euler_model = models.Euler(wave_speed_estimate_type = 'Davis')
# # Define left and right states for the test Riemann problem
# VL = np.array([1., 0., 1e+3])
# VR = np.array([1, 0., 1e-2])

# # Point of discontinuity 
# x0 = 0.5
# T_final = 0.012

# UL = euler_model.prim2cons(VL)
# UR = euler_model.prim2cons(VR)
# # Define function for initial conditions
# def u0_euler(x):
#     U0 = np.zeros((len(x), 3))
#     mask = x <= x0
#     U0[mask, :] = UL
#     mask = x > x0
#     U0[mask, :] = UR

#     return U0


# t1 = time.time()
# U, xs = solver(euler_model, nx = 200, nt=400, T = T_final, u0_fun = u0_euler,
#                 xl = 0., xr = + 1., bc_type = 'transparent')
# t2 = time.time()
# print ('Solution computed in {0:5.2e} seconds'.format(t2 - t1))

# plt.plot(xs, U[:, 0], 'o-')
# plt.grid(True)
# plt.show()

