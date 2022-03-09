#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
# Add your code

# Define function for initial conditions
def u0_2(x) :
    u0 = np.zeros((len(x), 1))
    u0[:, 0] = 1 + 0.5 * np.sin(np.pi*x)
    return u0

# Draw graph
def drawGraph(nx, nt, T) :
    fig = plt.figure(figsize = (16, 9))
    ax = fig.add_subplot(111)

    # Create model for burgers equation
    burgers_model = models.Burgers()

    # Compute numerical solution
    u, xs, dx, dt = solver(burgers_model, nx = nx, nt = nt, T = T, u0_fun = u0_2, xl = -1, xr = 1)
    
    # Plot initial condition
    xs0 = np.linspace(-1, 1, 100)
    ax.plot(xs0, u0_2(xs0), linestyle = '-', marker = 'None', color = 'black', label = 'Exact (T = 0)')
    
    # Plot numerical solution
    ax.plot(xs, u[:, 0], linestyle = '--', marker = 'x', markersize = 5, color = 'red', label = 'Numerical')
    
    
    # Configure plot
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
    
    plt.xlabel('x', fontsize = 18)
    plt.ylabel('u', fontsize = 18)

    plt.legend(loc = 'best', fontsize = 18)
    plt.title('dx = ' + str(dx) + ', dt = ' + str(dt) + ', T = ' + str(T), fontsize = 18)
    os.chdir(os.getcwd() + "\outData")
    plt.savefig('result dx = ' + str(dx) + ' dt = ' + str(dt) + ' T = ' + str(T) + '.jpg')
    os.chdir(os.getcwd().replace("\outData", ""))
    plt.show()
    
def drawAnimation(nx, nt) :
    fig = plt.figure(figsize = (16, 9))
    ax = fig.add_subplot(111)
    camera = Camera(fig)

    # Create model for burgers equation
    burgers_model = models.Burgers()

    # Compute initial condition
    xs0 = np.linspace(-1, 1, 100)
    ys0 = u0_2(xs0)
    
    # Create animaion
    tanim = np.linspace(0, 1.1, int(nt/5))
    for i in range(len(tanim)) : 
        u, xs, dx, dt = solver(burgers_model, nx = nx, nt = nt, T = tanim[i], u0_fun = u0_2, xl = -1, xr = 1)
        ax.plot(xs, u[:, 0], linestyle = '--', marker = 'x', color = 'red', markersize = 5)
        ax.plot(xs0, ys0, linestyle = '-', marker = 'None', color = 'black')
        plt.legend(['Numerical, T = ' + str(tanim[i]) + ' sec', 'Exact, T = 0 sec'], loc = 'best', fontsize = 18)
        camera.snap()
        plt.xlabel('x', fontsize = 18)
        plt.ylabel('u', fontsize = 18)
        plt.title('Animation ', fontsize = 18)
        
    # Configure plot
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
    
    animation = camera.animate()
    os.chdir(os.getcwd() + "\outData")
    animation.save('animation dx = ' + str(dx) + '.gif', fps = 20)
    os.chdir(os.getcwd().replace("\outData", ""))
    


# Configure input parameters (T, nx, nt)
T = 1.1
nx = 100
nt = 200

#drawGraph(nx, nt, 0)
#drawGraph(nx, nt, 0.5)
#drawGraph(nx, nt, 1.1)
#drawGraph(2 * nx, 4 * nt, 0)
#drawGraph(2 * nx, 4 * nt, 0.5)
#drawGraph(2 * nx, 4 * nt, 1.1)
drawAnimation(nx, nt)
#drawAnimation(2 * nx, 4 * nt)


# In[ ]:




