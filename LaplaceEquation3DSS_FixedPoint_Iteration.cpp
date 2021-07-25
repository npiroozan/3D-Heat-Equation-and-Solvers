//************************************

// Solution of the 3D Steady State Laplace Equation with the finite difference method - Fixed Point Iteration
// Code Developed by N. Piroozan

//************************************

#include <iostream>
#include <fstream>
#include <stdio.h>
using namespace std;

int main ()
{
    const int size = 50; // Grid Size (number of grid points)
    const int niter = 10000;  // Number of Iterations

    const double nx=50.0;     // Number of steps in the x direction
    const double ny=50.0;     // Number of steps in the y direction
    const double nz=50.0;     // Number of steps in the z direction

    const double dx=(20)/(nx-1);      // Width of space step (x)
    const double dy=(20)/(ny-1);      // Width of space step (y)
    const double dz=(20)/(nz-1);      // Width of space step (z)

    const double alpha1=dy * dy;
    const double alpha2=dx * dx;
    const double alpha3=dz * dz;

    const double beta1=alpha1*alpha3;
    const double beta2=alpha2*alpha3;
    const double beta3=alpha1*alpha2;
    const double beta4=2 * ((beta1) + (beta2) + (beta3));

    const int Tb=300;              // Initial temperature for the top edge

    double Temp[size][size][size]; 

        std::ofstream outfile ("Laplace3D_SS_Explicit.dat"); // Open an output file stream


    for (int k = 0; k < size; k++)  // Set up the boundary conditions
    {
        for (int j = 0; j < size; j++)
        {
            for (int i = 0; i < size; i++)
            {
                
                Temp[i][j][k]=0;
                Temp[i][j][0]=Tb;
            }
        }
    }


    for (int iter = 1; iter <= niter; iter++)
    {
        for (int i=1; i < (size-1); i++)    //x-direction
        {
            for (int j=1; j < (size-1); j++)
            {
                for (int k=1; k < (size-1); k++)
                {
                    Temp[i][j][k] = (beta1*(Temp[i+1][j][k] + Temp[i-1][j][k]) + beta2*(Temp[i][j+1][k] + Temp[i][j-1][k]) + beta3*(Temp[i][j][k+1] + Temp[i][j][k-1]))/(beta4);
                }
            }
        }
    }



    for (int k = 67; k < 68; k++)
    {
        for (int j = 0; j < size; j++)
        {
            for (int i = 0; i < size; i++)
            {
                outfile << Temp[i][j][2] << " ";
            }

            outfile << endl;
        }
        outfile << endl;
        
    }


    outfile.close ();
    return (0);
}