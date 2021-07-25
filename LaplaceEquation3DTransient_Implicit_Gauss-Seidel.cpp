//************************************

// Solution of the 3D Transient Laplace Equation with the finite difference method - Implicit Solution (Backwards Euler)
// Gauss-Seidel Solver
// Code Developed by N. Piroozan

//************************************

#include <iostream>
#include <fstream>
#include <stdio.h>
using namespace std;

int main ()
{
    const int size = 7; // Grid Size (number of grid points)
    const int nn = size*size*size;
    const int l = 30;   // Number of iterations for Gauss-Seidel

    const int nx=7;     // Number of steps in the x direction
    const int ny=7;     // Number of steps in the y direction
    const int nz=7;     // Number of steps in the z direction
    const int nt=7;     // Number of steps in the time domain

    const double dx=(2.0)/(nx-1);     // Size of each differential element in the x direction
    const double dy=(2.0)/(ny-1);     // Size of each differential element in the y direction
    const double dz=(2.0)/(nz-1);     // Size of each differential element in the z direction
    const double dt=(2.0)/(nt-1);     // Size of each Timestep

    const double Ti=300.0;          // Initial Temperature at t=0 for all x,y, and z
    const double kb=.003;           // Heat Conduction Coefficient (W/m*K)

    const double alpha0=kb*dt;      // Portion of calculation for Coefficients of the Global Stiffness Matrix
    const double alpha1=1/(dx*dx);  // Portion of calculation for Coefficients of the Global Stiffness Matrix
    const double alpha2=1/(dy*dy);  // Portion of calculation for Coefficients of the Global Stiffness Matrix
    const double alpha3=1/(dz*dz);  // Portion of calculation for Coefficients of the Global Stiffness Matrix

    double K[nn][nn]={0};           // Global Stiffness Matrix
    double F[nn]={0};               // Force Vector
    double T[nn]={0};               // Temperature Distribution
    double R[nn]={0};               // Residual

    std::ofstream outfile ("Laplace3D_Transient_Implicit_Gauss-Seidel.dat"); // Open an output file stream

    // Create the Global Stiffness Matrix

    for (int k=2; k <= (size-1); k++)
    {
        for (int j=2; j <= (size-1); j++)
        {
            for (int i=2; i <= (size-1); i++)
            {
                int nodenum = (i-1) + (j-1)*size + (k-1)*(size*size);

                K[nodenum][nodenum] = 1 + ((2*alpha0)*(alpha1 + alpha2 + alpha3));
                K[nodenum][nodenum+1] = -alpha0*alpha1;
                K[nodenum][nodenum-1] = -alpha0*alpha1;
                K[nodenum][nodenum+size] = -alpha0*alpha2;
                K[nodenum][nodenum-size] = -alpha0*alpha2;
                K[nodenum][nodenum+(size*size)] = -alpha0*alpha3;
                K[nodenum][nodenum-(size*size)] = -alpha0*alpha3;

            }
        }
    }

    // Set the Boundary Conditions

    // Set the Force Vector to reflect T=0 at all the surfaces of the system
    // Set the Global Stiffness Matrix Diagonal elements to value 1

    for (int k=1; k <= (size); k++)
    {
        for (int j=1; j <= (size); j++)
        {
            for (int i=1; i <= (size); i++)
            {
                int nodenum = (i-1) + (j-1)*size + (k-1)*(size*size);
                F[nodenum]=Ti;

                if (i==1)
                {
                    K[nodenum][nodenum]=1.0;
                    F[nodenum]=0;
                }

                else if(i==size)
                {
                    K[nodenum][nodenum]=1.0;
                    F[nodenum]=0;
                }

            }
        } 
    }

    for (int k=1; k <= (size); k++)
    {
        for (int j=1; j <= (size); j++)
        {
            for (int i=1; i <= (size); i++)
            {
                int nodenum = (i-1) + (j-1)*size + (k-1)*(size*size);

                if (j==1)
                {
                    K[nodenum][nodenum]=1.0;
                    F[nodenum]=0;
                }

                else if(j==size)
                {
                    K[nodenum][nodenum]=1.0;
                    F[nodenum]=0;
                }

            }
        } 
    }

    for (int k=1; k <= (size); k++)
    {
        for (int j=1; j <= (size); j++)
        {
            for (int i=1; i <= (size); i++)
            {
                int nodenum = (i-1) + (j-1)*size + (k-1)*(size*size);

                if (k==1)
                {
                    K[nodenum][nodenum]=1.0;
                    F[nodenum]=0;
                }

                else if(k==size)
                {
                    K[nodenum][nodenum]=1.0;
                    F[nodenum]=0;
                }


            }
        } 
    }


    // Solve for the Temperatures at throughout the system
    
    // Account for time here;  Temperatures calculated in each time step are used as the new Force Vector (F) in the next timestep
    // Global Stiffness Matrix (K) and the Force Vector (F)

 for (int p = 0; p < 7; p++)
    {

    //**********************************//
    // Gauss-Seidel Iterative Method to Solve for a System of Linear Equations
    for (int m = 0; m < l; m++)         // Number of Iterations
    {
        for (int i = 0; i < nn; i++)
        {
            R[i] = (F[i]/K[i][i]);            // Initialze the residual expression

            for (int j = 0; j < nn; j++)
            {
                if (i == j)
                {
                    continue;
                }

                else if (i != j)
                {
                    R[i] = R[i] - ((K[i][j] / K[i][i]) * T[j]);       // Set up the expressions for the residuals based on initialized solution vector
                    T[i] = R[i];                                      // Solution to the system of equations
                }
            }
        }
    }
    //**********************************//
        for (int i = 0; i < nn; i++)
        {
            F[i]=T[i];          // Solve for the updated solution vector and continue for each successive iteration
        }
    }

    //**********************************//
 

    // Output the Data File containing the Temperature Distribution at any Desired Timestep

    for (int i=0; i<nn; i++)
    {
        outfile << T[i] << " ";
        outfile << endl;
    }
    outfile << endl;

    outfile.close(); 
    return (0);
}
