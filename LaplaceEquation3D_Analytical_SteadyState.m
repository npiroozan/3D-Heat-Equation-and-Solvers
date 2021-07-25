% Solving 3D Steady State Laplace Equation from the Analytical Derivation
% Code written by N. Piroozan
clear;
clc;

%%
% Specify Parameters
nx=4;                         %Number of steps in the x direction
ny=4;                         %Number of steps in the y direction
nz=4;                         %Number of steps in the z direction

dx=2/(nx-1);                   %Size of each differential element in the x direction
dy=2/(ny-1);                   %Size of each differential element in the y direction
dz=2/(nz-1);                   %Size of each differential element in the z direction

x=0:dx:2;                      %Range of x(0,2) and specifying the grid points
y=0:dy:2;                      %Range of y(0,2) and specifying the grid points
z=0:dz:2;                      %Range of z(0,2) and specifying the grid points

W=2;                           %Size of the box in the x direction
L=2;                           %Size of the box in the y direction
H=2;                           %Size of the box in the z direction

M=130;                         %Number of series expansions for index i
N=130;                         %Number of series expansions for index j

Tb=300;                        %Temperature at z=L for all x and y

%%
% Analytical Solution

% Preallocating matrix T(i,j,k)
T = zeros(nx,ny,nz);

for i=1:nx
    for j=1:ny
        for k=1:nz
            for m = 1:M
                for n = 1:N
     
       kmn = sqrt(((m*pi)/W)^2+((n*pi)/L)^2);
 
        if (rem(m,2)~=0 && rem(n,2)~=0)
            Amn = (16*Tb)/(m*n*(pi^2))*(1/(sinh(kmn*H)));
        else
            Amn = 0;
        end
        
        T(i,j,k) = T(i,j,k) + Amn*sin(((m*pi)/W)*x(i))*sin(((n*pi)/L)*y(j))*sinh(kmn*z(k));

                end
            end
        end
    end
end

%%
% Postprocessing - Plot Results
Tp=T(:,:,3);

surf(x,y,Tp,'EdgeColor','none');
shading interp
title({'3D Steady State Laplace Equation - Analytical Solution'})
xlabel('X-axis (x) \rightarrow')
ylabel('{\leftarrow} Y-axis (y)')
zlabel('Temperature Profile (T) \rightarrow')            
            
