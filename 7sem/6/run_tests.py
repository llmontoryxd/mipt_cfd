#!/usr/bin/env python
# coding: utf-8

# In[2]:


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
import os
from celluloid import Camera

################################
# Test for advection equation
################################

# Create model for advection equation
#adv_model = models.Advection()

# Define function for initial conditions
#def u0_1(x):
#    u0 = np.zeros((len(x), 1))
#    u0[:, 0] = np.exp(-(x - 0.5)**2 / 0.01)
#    return u0

# Plot initial condition
#xs = np.linspace(0, 1, 50)
#plt.plot(xs, u0_1(xs), label = 'exact')

# Compute numerical solution
#t1 = time.time()
#u, xs = solver(adv_model, nx = 50, nt=100, T = 1, u0_fun = u0_1)
#t2 = time.time()
#print ('Solution computed in {0:5.2e} seconds'.format(t2 - t1))
# Plot numerical solution
#plt.plot(xs, u[:, 0], 'o', label = 'Numerical')
#plt.grid(True)
#plt.legend()
#plt.show()

################################
# Test for Burgers equation
################################

# Define function for initial conditions
#def u0_2(x) :
#    u0 = np.zeros((len(x), 1))
#    u0[:, 0] = 1 + 0.5 * np.sin(np.pi*x)
#    return u0

# Draw graph
#def drawGraph(nx, nt, T) :
#    fig = plt.figure(figsize = (16, 9))
#    ax = fig.add_subplot(111)

    # Create model for burgers equation
#    burgers_model = models.Burgers()

    # Compute numerical solution
#    u, xs, dx, dt = solver(burgers_model, nx = nx, nt = nt, T = T, u0_fun = u0_2, xl = -1, xr = 1)
    
    # Plot initial condition
#    xs0 = np.linspace(-1, 1, 100)
#    ax.plot(xs0, u0_2(xs0), linestyle = '-', marker = 'None', color = 'black', label = 'Exact (T = 0)')
    
    # Plot numerical solution
#    ax.plot(xs, u[:, 0], linestyle = '--', marker = 'x', markersize = 5, color = 'red', label = 'Numerical')
    
    
    # Configure plot
#    xax = ax.xaxis
#    xlocs = xax.get_ticklocs()
#    xlabels = xax.get_ticklabels()
#    xlines = xax.get_ticklines()

#    xax.grid(True)

#    for label in xlabels :
#        label.set_color('black')
#        label.set_rotation(0)
#        label.set_fontsize(18)
    
#    yax = ax.yaxis
#    ylocs = yax.get_ticklocs()
#    ylabels = yax.get_ticklabels()
#    ylines = yax.get_ticklines()
#    yax.grid(True)

#    for label in ylabels :
#        label.set_color('black')
#        label.set_rotation(0)
#        label.set_fontsize(18)
    
#    plt.xlabel('x', fontsize = 18)
#    plt.ylabel('u', fontsize = 18)

#    plt.legend(loc = 'best', fontsize = 18)
#    plt.title('dx = ' + str(dx) + ', dt = ' + str(dt) + ', T = ' + str(T), fontsize = 18)
#    os.chdir(os.getcwd() + "\outData")
#    plt.savefig('result dx = ' + str(dx) + ' dt = ' + str(dt) + ' T = ' + str(T) + '.jpg')
#    os.chdir(os.getcwd().replace("\outData", ""))
#    plt.show()
    
#def drawAnimation(nx, nt) :
#    fig = plt.figure(figsize = (16, 9))
#    ax = fig.add_subplot(111)
#    camera = Camera(fig)

    # Create model for burgers equation
#    burgers_model = models.Burgers()

    # Compute initial condition
#    xs0 = np.linspace(-1, 1, 100)
#    ys0 = u0_2(xs0)
    
    # Create animaion
#    tanim = np.linspace(0, 1.1, int(nt/5))
#    for i in range(len(tanim)) : 
#        u, xs, dx, dt = solver(burgers_model, nx = nx, nt = nt, T = tanim[i], u0_fun = u0_2, xl = -1, xr = 1)
#        ax.plot(xs, u[:, 0], linestyle = '--', marker = 'x', color = 'red', markersize = 5)
#        ax.plot(xs0, ys0, linestyle = '-', marker = 'None', color = 'black')
#        plt.legend(['Numerical, T = ' + str(tanim[i]) + ' sec', 'Exact, T = 0 sec'], loc = 'best', fontsize = 18)
#        camera.snap()
#        plt.xlabel('x', fontsize = 18)
#        plt.ylabel('u', fontsize = 18)
#        plt.title('Animation ', fontsize = 18)
        
    # Configure plot
#    xax = ax.xaxis
#    xlocs = xax.get_ticklocs()
#    xlabels = xax.get_ticklabels()
#    xlines = xax.get_ticklines()

#    xax.grid(True)

#    for label in xlabels :
#        label.set_color('black')
#        label.set_rotation(0)
#        label.set_fontsize(18)
    
#    yax = ax.yaxis
#    ylocs = yax.get_ticklocs()
#    ylabels = yax.get_ticklabels()
#    ylines = yax.get_ticklines()
#    yax.grid(True)

#    for label in ylabels :
#        label.set_color('black')
#        label.set_rotation(0)
#        label.set_fontsize(18)
    
#    animation = camera.animate()
#    os.chdir(os.getcwd() + "\outData")
#    animation.save('animation dx = ' + str(dx) + '.gif', fps = 20)
#    os.chdir(os.getcwd().replace("\outData", ""))
    


# Configure input parameters (T, nx, nt)
#T = 1.1
#nx = 100
#nt = 200

#drawGraph(nx, nt, 0)
#drawGraph(nx, nt, 0.5)
#drawGraph(nx, nt, 1.1)
#drawGraph(2 * nx, 4 * nt, 0)
#drawGraph(2 * nx, 4 * nt, 0.5)
#drawGraph(2 * nx, 4 * nt, 1.1)
#drawAnimation(nx, nt)
#drawAnimation(2 * nx, 4 * nt)


###############################
# Tests for Euler equations
###############################

def drawGraphEuler(rs_type, VL, VR, nx, nt, x0, xl, xr, T, bc_type, name) : 
    # define models for Davis and Einfeldt evaluations
    euler_modelDavis = models.Euler(rs_type = rs_type, wave_speed_estimate_type = 'Davis')
    euler_modelEinfeldt = models.Euler(rs_type = rs_type, wave_speed_estimate_type = 'Einfeldt')
    
    # UL and UR - vectors of conservative variables from VL and VR
    UL = euler_modelDavis.prim2cons(VL)
    UR = euler_modelDavis.prim2cons(VR)
    
    # get solution
    UDavis, xsDavis = solver(model = euler_modelDavis, nx = nx, nt = nt, u0_fun = lambda x: u0_euler(x, x0, UL, UR), xl = xl, xr = xr, T = T, bc_type = bc_type)
    UEinfeldt, xsEinfeldt = solver(model = euler_modelEinfeldt, nx = nx, nt = nt, u0_fun = lambda x: u0_euler(x, x0, UL, UR), xl = xl, xr = xr, T = T, bc_type = bc_type)
    
    # VDavis and VEinfeldt - vectors of primitive variables from UDavis and UEinfeldt
    VDavis = np.array([euler_modelDavis.cons2prim(x) for x in UDavis])
    VEinfeldt = np.array([euler_modelDavis.cons2prim(x) for x in UEinfeldt])
    
    # draw graphs
    fig = plt.figure(figsize = (16, 9))
    
    ax = fig.add_subplot(221)
    ax.plot(xsDavis, VDavis[:, 0], linestyle = '--', color = 'black', marker = 'x', markersize = 3, label = 'Davis')
    ax.plot(xsEinfeldt, VEinfeldt[:, 0], linestyle = '--', color = 'red', marker = 'o', markersize = 3, label = 'Einfeldt') 
    xax = ax.xaxis
    xlocs = xax.get_ticklocs()
    xlabels = xax.get_ticklabels()
    xlines = xax.get_ticklines()
    xax.grid(True)
    for label in xlabels :
        label.set_color('black')
        label.set_rotation(0)
        label.set_fontsize(18) 
    yax = ax.yaxis
    ylocs = yax.get_ticklocs()
    ylabels = yax.get_ticklabels()
    ylines = yax.get_ticklines()
    yax.grid(True)
    for label in ylabels :
        label.set_color('black')
        label.set_rotation(0)
        label.set_fontsize(18)
    ax.set_xlabel('x', fontsize = 18)
    ax.set_ylabel('Density', fontsize = 18)
    ax.legend(loc = 'best', fontsize = 18)
    
    ax = fig.add_subplot(222)
    ax.plot(xsDavis, VDavis[:, 1], linestyle = '--', color = 'black', marker = 'x', markersize = 3, label = 'Davis')
    ax.plot(xsEinfeldt, VEinfeldt[:, 1], linestyle = '--', color = 'red', marker = 'o', markersize = 3, label = 'Einfeldt') 
    xax = ax.xaxis
    xlocs = xax.get_ticklocs()
    xlabels = xax.get_ticklabels()
    xlines = xax.get_ticklines()
    xax.grid(True)
    for label in xlabels :
        label.set_color('black')
        label.set_rotation(0)
        label.set_fontsize(18) 
    yax = ax.yaxis
    ylocs = yax.get_ticklocs()
    ylabels = yax.get_ticklabels()
    ylines = yax.get_ticklines()
    yax.grid(True)
    for label in ylabels :
        label.set_color('black')
        label.set_rotation(0)
        label.set_fontsize(18)
    ax.set_xlabel('x', fontsize = 18)
    ax.set_ylabel('Velocity', fontsize = 18)
    ax.legend(loc = 'best', fontsize = 18)
    
    ax = fig.add_subplot(223)
    ax.plot(xsDavis, VDavis[:, 2], linestyle = '--', color = 'black', marker = 'x', markersize = 3, label = 'Davis')
    ax.plot(xsEinfeldt, VEinfeldt[:, 2], linestyle = '--', color = 'red', marker = 'o', markersize = 3, label = 'Einfeldt') 
    xax = ax.xaxis
    xlocs = xax.get_ticklocs()
    xlabels = xax.get_ticklabels()
    xlines = xax.get_ticklines()
    xax.grid(True)
    for label in xlabels :
        label.set_color('black')
        label.set_rotation(0)
        label.set_fontsize(18) 
    yax = ax.yaxis
    ylocs = yax.get_ticklocs()
    ylabels = yax.get_ticklabels()
    ylines = yax.get_ticklines()
    yax.grid(True)
    for label in ylabels :
        label.set_color('black')
        label.set_rotation(0)
        label.set_fontsize(18)
    ax.set_xlabel('x', fontsize = 18)
    ax.set_ylabel('Pressure', fontsize = 18)
    ax.legend(loc = 'best', fontsize = 18)
    
    
    os.chdir(os.getcwd() + "\outData")
    plt.savefig('Test ' + name + ', method = ' + rs_type + '.jpg')
    os.chdir(os.getcwd().replace("\outData", ""))
    
# Define function for initial conditions
def u0_euler(x, x0, UL, UR):
    U0 = np.zeros((len(x), 3))
    mask = x <= x0
    U0[mask, :] = UL
    mask = x > x0
    U0[mask, :] = UR

    return U0
    
# define nx and nt
nx = 400
nt = 1200

# define xl and xr
xl = 0.0
xr = 1.0

##################################
## Test 1
##################################
print("Test 1 starts")

x0 = 0.3
T_final = 0.2
VL = np.array([1, 0.75, 1])
VR = np.array([0.125, 0.0, 0.1])

drawGraphEuler('hll', VL, VR, nx, nt, x0, xl, xr, T_final, "transparent", "1")
drawGraphEuler('hllc', VL, VR, nx, nt, x0, xl, xr, T_final, "transparent", "1")

print("Test 1 completed!")

##################################
## Test 2
##################################
print("Test 2 starts")

x0 = 0.5
T_final = 0.15
VL = np.array([1., -2.0, 0.4])
VR = np.array([1., 2.0, 0.4])

drawGraphEuler('hll', VL, VR, nx, nt, x0, xl, xr, T_final, "transparent", "2")
drawGraphEuler('hllc', VL, VR, nx, nt, x0, xl, xr, T_final, "transparent", "2")

print("Test 2 completed!")

##################################
## Test 3
##################################
print("Test 3 starts")

x0 = 0.5
T_final = 0.012
VL = np.array([1., 0., 1000.])
VR = np.array([1., 0., 0.01])

drawGraphEuler('hll', VL, VR, nx, nt, x0, xl, xr, T_final, "transparent", "3")
drawGraphEuler('hllc', VL, VR, nx, nt, x0, xl, xr, T_final, "transparent", "3")

print("Test 3 completed!")

##################################
## Test 4
##################################
print("Test 4 starts")

x0 = 0.4
T_final = 0.035
VL = np.array([5.99924, 19.5975, 460.894])
VR = np.array([5.99242, -6.19633, 46.0950])

drawGraphEuler('hll', VL, VR, nx, nt, x0, xl, xr, T_final, "transparent", "4")
drawGraphEuler('hllc', VL, VR, nx, nt, x0, xl, xr, T_final, "transparent", "4")

print("Test 4 completed!")

##################################
## Test 5
##################################
print("Test 5 starts")

x0 = 0.8
T_final = 0.012
VL = np.array([1., -19.59745, 1000.])
VR = np.array([1., -19.59745, 0.01])

drawGraphEuler('hll', VL, VR, nx, nt, x0, xl, xr, T_final, "transparent", "5")
drawGraphEuler('hllc', VL, VR, nx, nt, x0, xl, xr, T_final, "transparent", "5")

print("Test 5 completed!")


# In[ ]:




