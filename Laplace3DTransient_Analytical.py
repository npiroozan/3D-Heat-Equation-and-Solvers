##****************************************

# Analytical Solution of the 3D Laplace Equation
# Code Developed in Python by N. Piroozan

##****************************************

import matplotlib.pyplot as plt
import numpy as np
import math

nx = 4                             # Number of steps in the x direction
ny = 4                             # Number of steps in the y direction
nz = 4                             # Number of steps in the z direction
nt = 4                             # Number of steps in the time domain

dx = (2)/(nx-1)                     # Width of space step (x)
dy = (2)/(ny-1)                     # Width of space step (y)
dz = (2)/(nz-1)                     # Width of space step (z)
dt = (2)/(nt-1)                     # Size of each Timestep (t)

x = np.arange(0,2+.1,dx)            # Range of x(0,2) and specifying grid points
y = np.arange(0,2+.1,dy)            # Range of y(0,2) and specifying grid points
z = np.arange(0,2+.1,dz)            # Range of z(0,2) and specifying grid points
t = np.arange(0,2+.1,dt)            # Range of time domain

W = 2                               # Total length of the x-axis
L = 2                               # Total length of the y-axis
H = 2                               # TOtal length of the z-axis

M = 60                             # Number of taylor series expansions for index i
N = 60                             # Number of taylor series expansions for index j
P = 60                             # Number of taylor series expansions for index k

Ti = 300                            # Constant Temperature for Z=H at all x and y
kb = 0.003                          # Heat Conduction Coefficient in W/(m*K)

Temp = np.zeros((nx,ny,nz,nt))      # Preallocating matrix T(i,j,k,l)

for l in range(nt):
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for m in range(1,M):
                    for n in range(1,N):
                        for p in range(1,P):

                            Amnl = ((((2*m-1)*math.pi)/W)**2) + ((((2*n-1)*math.pi)/L)**2) + ((((2*p-1)*math.pi)/H)**2)
                            mum = ((2*m-1)*math.pi)/W
                            vun = ((2*n-1)*math.pi)/L
                            kal = ((2*p-1)*math.pi)/H

                            Temp(i,j,k,l) = Temp(i,j,k,l) + ((64*(Ti))/((math.pi)**3))*(((math.sin(mum*x(i)))*(math.sin(vun*y(j)))*((math.sin(kal*z(k))))*math.exp(-Amnl*kb*t(l)))/((2*m-1)*(2*n-1)*(2*p-1)))


np.savetxt("Laplace3DTransient_Analytical.txt", Temp, fmt="%s")
