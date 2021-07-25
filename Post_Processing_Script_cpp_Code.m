%% Post-Processing Script for C++ Data File
%  Explicit Solution
clear;
clc;


%% Load Data file, called Laplace3D_Transient
load Laplace3D_Transient.dat;


%% Specify the Parameters
%  Same as the parameters used for the C++ Simulation

% Specify Parameters
nx=20;                         %Number of steps in the x direction
ny=20;                         %Number of steps in the y direction 
nz=20;                         %Number of steps in the z direction
nt=20;                         %Number of steps in the time domain

dx=2/(nx-1);                   %Size of each differential element in the x direction
dy=2/(ny-1);                   %Size of each differential element in the y direction
dz=2/(nz-1);                   %Size of each differential element in the z direction
dt=2/(nt-1);                   %Size of each Timestep

x=0:dx:2;                      %Range of x(0,2) and specifying the grid points
y=0:dy:2;                      %Range of y(0,2) and specifying the grid points
z=0:dz:2;                      %Range of z(0,2) and specifying the grid points
t=0:dt:2;                      %Range of time domain


%% Plot the Results for the data set
%  In this example, we are plotting for all x and y and at a particular z
%  and t

surf(x,y,Laplace3D_Transient,'EdgeColor','none');
shading interp
title({'3D Transient Heat Conduction Equation - Explicit Solution from C++ Code'})
xlabel('X-axis (x) \rightarrow')
ylabel('{\leftarrow} Y-axis (y)')
zlabel('Temperature Profile (T) \rightarrow')            
