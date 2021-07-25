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

dx = (2)/(nx-1)                     # Width of space step (x)
dy = (2)/(ny-1)                     # Width of space step (y)
dz = (2)/(nz-1)                     # Width of space step (z)

x = np.arange(0,2+.1,dx)
y = np.arange(0,2+.1,dy)
z = np.arange(0,2+.1,dz)

W = 2                               # Total length of the x-axis
L = 2                               # Total length of the y-axis
H = 2                               # TOtal length of the z-axis

M = 100                              # Number of taylor series expansions for x
N = 100                              # Number of taylor series expansions for y

Tb = 300                            # Constant Temperature for Z=H at all x and y

Temp = np.zeros((nx,ny,nz))         # Preallocating matrix T(i,j,k)


# Solution of the Analytical Formulation

for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            for m in range(1,M):
                for n in range(1,N):

                    kmn = math.sqrt(((m*math.pi)/W)**2 + ((n*math.pi)/L)**2)

                    if (m&1==1 and n&1==1):

                        Amn = (16*Tb)/(m*n*(math.pi)**2)*(1/(math.sinh(kmn*H)))

                    else:

                        Amn = 0

                    Temp[i][j][k] = Temp[i][j][k] + Amn*math.sin(((m*math.pi)/W)*x[i])*math.sin(((n*math.pi)/L)*y[j])*math.sinh(kmn*z[k])


np.savetxt("Laplace3DSS_Analytical.txt", Temp[:][:][2], fmt="%s")

