//************************************

// Solver for a System of Simultaneous Linear Equations using Gauss-Jordan Elimination Method
// Code Developed by N. Piroozan

//************************************

#include <iostream>
#include <fstream>
#include <stdio.h>
using namespace std;

int main()
{
    int nn;

    std::ofstream outfile("Gauss-Jordan_Solver_Results.dat"); // Open an output file stream

    cout<<"\n Enter the number of simultaneous linear equations you wish to solve for: ";
    cin >> nn;

    // If the number of simultaneous linear equations is nn, then the size of the augmented matrix, Aug, will be nn*(nn+1)
    double Aug[nn][nn+1];

    // For nn equations, we will have nn unknowns which will be stored within the array called Result (Res), of size nn
    double res[nn];

    // In the script "LaplaceEquation3DSS_Implicit.cpp", the Augmented matrix is determined within the code.
    // For this standalone code, please enter the elements of the Augmented matrix manually
    cout << "\n Enter the elements of the Augmented Matrix: ";

    for (int i = 0; i < nn; i++)
     {
        for (int j = 0; j < (nn+1); j++)
        {
            cin >> Aug[i][j];
        }
     }


    //double Aug[3][4] = { {1,0,0,300},{2,-10,-10,100},{3,15,25,200} };

    printf("Original Augmented Matrix: \n");
    for (int p = 0; p < nn; p++)
    {
        for (int q = 0; q < nn + 1; q++)
        {
            printf("%10.3f", Aug[p][q]);
        }
        printf("\n");
    }

    printf("\n\n Swapping Rows: \n");

    // Solver for the System of Linear Equations
    for (int i = 0; i < nn-1; i++)
    {
        for (int j = i + 1; j < nn; j++)
        {
            if (abs(Aug[i][i]) < abs(Aug[j][i]))
            {
                for (int k = 0; k < nn + 1; k++)
                {
                    // Swap Aug(i,k) and Aug(j,k)
                    Aug[i][k] = Aug[i][k] + Aug[j][k];
                    Aug[j][k] = Aug[i][k] - Aug[j][k];
                    Aug[i][k] = Aug[i][k] - Aug[j][k];
                }
            }
        }
        for (int p = 0; p < nn; p++)
        {
            for (int q = 0; q < nn + 1; q++)
            {
                printf("%10.3f", Aug[p][q]);
            }
            printf("\n");
        }
        printf("\n");
    }

    printf("\n Gaussian Elimination of the Lower Triangle: \n");
    // Perform Gaussian Elimination of the lower triangle in the Augmented Matrix
    for (int i = 0; i < nn - 1; i++)
    { 
        // Normalizing and Pivoting
        float f1 = 1. / Aug[i][i];
        for (int k = 0; k < nn + 1; k++)
        {
            Aug[i][k] = f1 * Aug[i][k];
        }   
        for (int j = i + 1; j < nn; j++)
        {
            float f2 = Aug[j][i];
            for (int k = 0; k < nn + 1; k++)
            {
                Aug[j][k] = Aug[j][k] - f2 * Aug[i][k];
            }
        }
    }

    for (int i = 0; i < nn; i++)
    {
        for (int j = 0; j < nn + 1; j++)
        {
            printf("%10.3f", Aug[i][j]);
        }
        printf("\n");
    }

    printf("\n\n Gaussian Elimination of the Upper Triangle: \n");
    // Perform Gaussian Elimination of the upper triangle in the Augmented Matrix
    for (int i = nn-1; i >= 0; i--)
    {
        // Normalizing and Pivoting
         float f3 = 1. / Aug[nn-1][nn-1];
         for (int k = nn-1; k < nn + 1; k++)
         {
             Aug[nn-1][k] = f3 * Aug[nn-1][k];
         }
        for (int p = 0; p < nn; p++)
        {
            for (int q = 0; q < nn + 1; q++)
            {
                printf("%10.3f", Aug[p][q]);
            }
            printf("\n");
        }
        for (int j = i-1; j >= 0; j--)
        {
            float f4 = Aug[j][i];
            for (int k = i; k < nn + 1; k++)
            {
                Aug[j][k] = Aug[j][k] - f4 * Aug[i][k];
            }
        }
        printf("\n");
    }

    printf("\n Results:");
    for (int i = 0; i < nn; i++)
    {
        printf("\n %10.3f", Aug[i][nn]);
    }
    printf("\n");

    for (int i = 0; i < nn; i++)
    {
        outfile << Aug[i][nn] << " ";
        outfile << endl;
    }
    outfile << endl;
    outfile.close();       
    return (0);
}
