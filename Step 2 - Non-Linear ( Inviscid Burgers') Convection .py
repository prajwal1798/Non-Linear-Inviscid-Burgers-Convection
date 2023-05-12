#!/usr/bin/env python
# coding: utf-8

# ## STEP 2 :       1-D Non-Linear (Inviscid Burgers') Convection

# Lets start by importing necessary libraries 

# In[3]:


import numpy as np
import matplotlib.pyplot as plt
from math import *


# Defining the domain for both Space & Time 

# In[4]:


xmax = 2                         # Total Domain lenth in (m)
nx   = 41                       # Number of Grid Points in space 'X'
dx   = xmax/(nx-1)               # Size of each grid in (m)
nt   = 25                        # Number of Grid Points in time 't'
dt   = 0.025                     # Time-step size in (s)
c    = 1                         # Wave Propagation Speed in (m/s)
x    = np.linspace(0,xmax,nx)    # Space Domain with grids in (m)


# Defining the Initial Conditions

# In[5]:


u = np.ones(nx)                      # Initialise the velocity array of ones
u[int(.5/dx):int((1/dx)+1)] = 2      # Implementing the square wave condition for velocity


# Plotting the Velocity field function U to understand its profile variation in Space

# In[6]:


plt.plot(x,u)
plt.title('Velocity Profile - Square Wave')
plt.xlabel('X (m)')
plt.ylabel('u (m/s)')
plt.grid()
plt.show


# Discretisation of Linear Convection equation & finding its max grid resolution to avoid solution blow-up  
# 
# $$\frac{\partial u}{\partial t} + u* \frac{\partial u}{\partial x} = 0$$ 
# 
# Forward Differencing in Time :
# $$\frac{\partial u}{\partial t} \approx \frac{u_i^{n+1}-u_i^{n}}{\Delta t}\rightarrow 1$$
# 
# Backward Differencing in Space :
# $$\frac{\partial u}{\partial x} \approx \frac{u_i^{n}-u_{i-1}^{n}}{\Delta x}\rightarrow 2$$
# 
# Therefore expressing Velocity explicitly, we have;
# $$u_i^{n+1} = u_i^{n} -u_i^{n} *\frac{\Delta t}{\Delta x}*(u_i^{n}-u_{i-1}^{n}) \rightarrow 3$$  

# Lets start by Intialising a new array of velocity $u_n$

# In[7]:


#Initialize a temporary array
un = np.ones(nx) 

# Time Loop
for n in range(nt):
    un = u.copy()
# Space Loop
    for i in range(1,nx):
        u[i] = un[i]-un[i]*dt/dx*(un[i]-un[i-1])


# Lets Plot the following

# In[8]:


plt.plot(x,u)
plt.title('Convected Profile Velocity')
plt.xlabel('X (m)')
plt.ylabel('u (m/s)')
plt.grid()
plt.show

