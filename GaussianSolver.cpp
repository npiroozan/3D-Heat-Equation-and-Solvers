//************************************

// Solver for a System of Simultaneous Linear Equations using Gaussian Elimination
// Code Developed by N. Piroozan

//************************************

#include <iostream>
#include <fstream>
#include <stdio.h>
using namespace std;

int main ()
{

    int nn;

    std::ofstream outfile ("Gaussian_Solver_Results.dat"); // Open an output file stream

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

   // Solver for the System of Linear Equations

    for(int i=0;i<nn;i++) 
    {                   
        for(int j=i+1;j<nn;j++)
        {
            if(abs(Aug[i][i]) < abs(Aug[j][i]))
            {
                for(int k=0;k<nn+1;k++)
                {
                    // Swap Aug(i,k) and Aug(j,k)
                    Aug[i][k]=Aug[i][k]+Aug[j][k];
                    Aug[j][k]=Aug[i][k]-Aug[j][k];
                    Aug[i][k]=Aug[i][k]-Aug[j][k];
                }
            }
      }
    }
   
     // Perform Gaussian Elimination of the lower triangle in the Augmented Matrix
    for(int i=0;i<nn-1;i++)
    {
        for(int j=i+1;j<nn;j++)
        {
            float f=Aug[j][i]/Aug[i][i];
            for(int k=0;k<nn+1;k++)
            {
              Aug[j][k]=Aug[j][k]-f*Aug[i][k];
      }
        }
    }

    // Perform backwards substitution to determine the Temperatures
    for(int i=nn-1;i>=0;i--)          
    {                     
        res[i]=Aug[i][nn];
                    
        for(int j=i+1;j<nn;j++)
        {
          if(i!=j)
          {
              res[i]=res[i]-Aug[i][j]*res[j];
    }          
  }
  res[i]=res[i]/Aug[i][i];  
    }

    for (int i=0; i<nn; i++)
    {
        outfile << res[i] << " ";
        outfile << endl;
    }
    outfile << endl;

    

     outfile.close(); 
     return (0);

}