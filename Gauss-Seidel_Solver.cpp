//************************************

// Solver for a System of Simultaneous Linear Equations using the Gauss-Seidel Iteration Method
// Code Developed by N. Piroozan

//************************************

#include <iostream>
#include <fstream>
#include <stdio.h>
#include<cmath>
using namespace std;

int main()
{
    int i,j;        // General Indicies to be used for System of Equation loops
    int m,l;        // General Indices used for Iteration Loop
    int nn;         // Number of simultaneous linear equations

    std::ofstream outfile("Gauss-Seidel_Solver.out"); // Open an output file stream

    // Enter the value for the number of linear equations (nn)
    cout<<"\n Enter the number of linear equations: \n";
    cin>>nn;

    double T[nn],G[nn][nn],F[nn],R[nn];         // Size of the vectors and matrix in the solver

    // T[nn] refers to the size of the solution array
    // G[nn][nn] refers to the size of the Global Stiffness Matrix
    // F[nn] refers to the size of the force vector
    // R[nn] refers to the size of the residual vector

    // Enter the number of iterations (l)
    cout<<"\n Enter the desired number of iterations: \n";
    cin>>l;

    // Enter the values of in the force vector (F[nn])
    cout<<"\n Enter the values in the Force vector: \n";

    for(i = 0; i < nn; i++)
    {
        cin>>F[i];
    }

    // Coefficients of the Global Stiffnes Matrix (a[nn][nn])
    cout<<"\n Enter the coefficients of the Global Stiffness Matrix: \n";

    for(i = 0; i < nn; i++) 
    {
        T[i] = 0;                   // Initialze the solution vector to 0
    
        for(j = 0; j < nn; j++) 
        {
            cin >> G[i][j];
        }
    }


    //**********************************//
    // Gauss-Seidel Iterative Method to Solve for a System of Linear Equations
    for (m = 0; m < l; m++)         // Number of Iterations
    {
        for (i = 0; i < nn; i++)
        {
            R[i] = (F[i]/G[i][i]);            // Initialze the residual expression

            for (j = 0; j < nn; j++)
            {
                if (i == j)
                {
                    continue;
                }

                else if (i != j)
                {
                    R[i] = R[i] - ((G[i][j] / G[i][i]) * T[j]);       // Set up the expressions for the residuals based on initialized solution vector
                    T[i] = R[i];                                      // Solution to the system of equations
                }
            }
        }
    }
    //**********************************//


    // Save the results in an output file
    for (int i = 0; i < nn; i++)
    {
        outfile << T[i] << " ";
        outfile << endl;
    }
 
    outfile.close();
    return (0);
}
