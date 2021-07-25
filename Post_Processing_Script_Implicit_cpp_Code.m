%% Post-Processing Script for C++ Data File
%  Explicit Solution
clear;
clc;


%% Load Data file, called Laplace3D_Transient
load Laplace3DSS.dat;


%% Specify the Parameters
%  Same as the parameters used for the C++ Simulation

% Specify Parameters
nx=7;                         %Number of steps in the x direction
ny=7;                         %Number of steps in the y direction 
nz=7;                         %Number of steps in the z direction
nt=7;                         %Number of steps in the time domain

dx=2/(nx-1);                   %Size of each differential element in the x direction
dy=2/(ny-1);                   %Size of each differential element in the y direction
dz=2/(nz-1);                   %Size of each differential element in the z direction
dt=2/(nt-1);                   %Size of each Timestep

x=0:dx:2;                      %Range of x(0,2) and specifying the grid points
y=0:dy:2;                      %Range of y(0,2) and specifying the grid points
z=0:dz:2;                      %Range of z(0,2) and specifying the grid points
t=0:dt:2;                      %Range of time domain


%% Convert Node Number, of the C++ Results, to Cartesian Coordinate System
% Please recall in the C++ code, node number was used to convert the 3D
% system into a new basis to allow for a 2D Global Stiffness Matrix and
% thereby be able to solve for a linear system of equations

for k=1:nz
    for j=1:ny
        for i=1:nx
        
        nodenum = i+(j-1)*nx+(k-1)*nx*ny        %Convert i,j,k into node number
        Tn(i,j,k) = Laplace3DSS(nodenum);
        
        end 
    end
end

% Select the desired cross section to plot (the particular z)
Tp = Tn (:,:,5);


%% Plot the Results for the data set
%  In this example, we are plotting for all x and y and at a particular z
%  and t

surf(x,y,Tp,'EdgeColor','none');
shading interp
title({'3D Transient Heat Conduction Equation - Implicit Solution from the C++ Code'})
xlabel('X-axis (x) \rightarrow')
ylabel('{\leftarrow} Y-axis (y)')
zlabel('Temperature Profile (T) \rightarrow')            
