//************************************

// Solution of the 3D Steady State Laplace Equation with the finite difference method - Explicit Method (Jacobi Solver)
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

    const double nx=7.0;     // Number of steps in the x direction
    const double ny=7.0;     // Number of steps in the y direction
    const double nz=7.0;     // Number of steps in the z direction

    const double dx=(2)/(nx-1);     // Width of space step (x)
    const double dy=(2)/(ny-1);     // Width of space step (y)
    const double dz=(2)/(nz-1);     // Width of space step (z)

    const double alpha1=dy * dy;
    const double alpha2=dx * dx;
    const double alpha3=dz * dz;

    const double beta1=alpha1*alpha3;
    const double beta2=alpha2*alpha3;
    const double beta3=alpha1*alpha2;

    const double Tb=300.0;

    double K[nn][nn]={0};
    double F[nn]={0};
    double T[nn]={0};
    double R[nn]={0};

    std::ofstream outfile ("Laplace3DSS_Explicit_Jacobi.dat"); // Open an output file stream


    for (int k=2; k <= (size-1); k++)
    {
        for (int j=2; j <= (size-1); j++)
        {
            for (int i=2; i <= (size-1); i++)
            {
                int nodenum = (i-1) + (j-1)*size + (k-1)*(size*size);

                K[nodenum][nodenum] = -2*(beta1 + beta2 + beta3);
                K[nodenum][nodenum+1] = beta1;
                K[nodenum][nodenum-1] = beta1;
                K[nodenum][nodenum+size] = beta2;
                K[nodenum][nodenum-size] = beta2;
                K[nodenum][nodenum+(size*size)] = beta3;
                K[nodenum][nodenum-(size*size)] = beta3;
     
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

                if (i==1)
                {
                    K[nodenum][nodenum]=1.0;
                }

                else if(i==size)
                {
                    K[nodenum][nodenum]=1.0;
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
                }

                else if(j==size)
                {
                    K[nodenum][nodenum]=1.0;
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
                }

                else if(k==size)
                {
                    K[nodenum][nodenum]=1.0;
                    F[nodenum]=Tb;
                }


            }
        } 
    }


    //**********************************//
    // Jacobi Iterative Method to Solve for a System of Linear Equations
    for (int m = 0; m < l; m++)         // Number of Iterations
    {
        for (int i = 0; i < nn; i++)
        {
            R[i] = F[i];            // Initialze the residual expression to be equal to the Force vector

            for (int j = 0; j < nn; j++)
            {
                if (i != j)
                {
                    R[i] = R[i] - K[i][j] * T[j];       // Set up the expressions for the residuals based on initialized solution vector
                }
            }
        }

        for (int i = 0; i < nn; i++)
        {
            T[i] = R[i] / K[i][i];          // Solve for the updated solution vector and continue for each successive iteration
        }
    }
    //**********************************//


    for (int i=0; i<nn; i++)
    {
        outfile << T[i] << " ";
        outfile << endl;
    }
    outfile << endl;

     outfile.close(); 
     return (0);
}
