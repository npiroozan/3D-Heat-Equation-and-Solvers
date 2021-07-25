//************************************

// Solution of the 3D Transient Laplace Equation with the finite difference method - Fixed Point Iteration
// Code Developed by N. Piroozan

//************************************

#include <iostream>
#include <fstream>
#include <stdio.h>
using namespace std;

int main ()
{
    const int size = 20; // Grid Size (number of grid points)
    const int niter = 10000;  // Number of Iterations

    const double nx=20;     // Number of steps in the x direction
    const double ny=20;     // Number of steps in the y direction
    const double nz=20;     // Number of steps in the z direction
    const double nt=20;     // Number of steps in the time domain

    const double dx=(2)/(nx-1);      // Width of space step (x)
    const double dy=(2)/(ny-1);      // Width of space step (y)
    const double dz=(2)/(nz-1);      // Width of space step (z)
    const double dt=(2)/(nt-1);       // Size of each Time step

    const double alpha1=dy * dy;
    const double alpha2=dx * dx;
    const double alpha3=dz * dz;

    const int Ti=300;              // Initial temperature for the top edge
    const double kb=0.003;          // Heat Conduction Coefficient (W/m*K)

    const double beta1 = 1 - 2*kb*dt*((1/alpha2)+(1/alpha1)+(1/alpha3));
    const double beta2=(kb*dt)/(alpha2);
    const double beta3=(kb*dt)/(alpha1);
    const double beta4=(kb*dt)/(alpha3);

    double Temp[size][size][size][size]; 

        std::ofstream outfile ("Laplace3D_Transient.dat"); // Open an output file stream


        for (int l = 0; l < size ; l++)
        {
            for (int k = 0; k < size; k++)
            {
                for (int j = 0; j < size; j++)
                {
                    for (int i = 0; i < size; i++)
                    {

                        Temp[i][j][k][l]=0;
                        Temp[i][j][k][0]=Ti;

                    }
                }
            }
        }

        for (int iter = 1; iter < niter; iter++)
        {
            for (int l = 0; l < (size-1); l++)
            {
                for (int i = 1; i < (size-1); i++)
                {
                    for (int j = 1; j < (size-1); j++)
                    {
                        for (int k = 1; k < (size-1); k++)
                        {
                            Temp[i][j][k][l+1] = (Temp[i][j][k][l])*beta1 + beta2*(Temp[i+1][j][k][l]+Temp[i-1][j][k][l]) + beta3*(Temp[i][j+1][k][l]+Temp[i][j-1][k][l]) + beta4*(Temp[i][j][k+1][l]+Temp[i][j][k-1][l]);
                        }
                    }
                }
            }
        }


        for (int l = 17; l < 18; l++)
        {
            for (int k = 17; k < 18; k++)
            {
                for (int j = 0; j < size; j++)
                {
                    for (int i = 0; i < size; i++)
                    {
                        outfile << Temp[i][j][17][17] << " ";
                    }

                    outfile << endl;
                }

                outfile << endl;
            }

            outfile << endl;
        }

        outfile.close();
        return (0);

}
