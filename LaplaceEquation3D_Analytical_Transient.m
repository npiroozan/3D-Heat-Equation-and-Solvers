% Solving 3D Transient Laplace Equation from the Analytical Derivation
% Code written by N. Piroozan
clear;
clc;

%%
% Specify Parameters
nx=4;                         %Number of steps in the x direction
ny=4;                         %Number of steps in the y direction 
nz=4;                         %Number of steps in the z direction
nt=4;                         %Number of steps in the time domain

dx=2/(nx-1);                   %Size of each differential element in the x direction
dy=2/(ny-1);                   %Size of each differential element in the y direction
dz=2/(nz-1);                   %Size of each differential element in the z direction
dt=2/(nt-1);                 %Size of each Timestep

x=0:dx:2;                      %Range of x(0,2) and specifying the grid points
y=0:dy:2;                      %Range of y(0,2) and specifying the grid points
z=0:dz:2;                      %Range of z(0,2) and specifying the grid points
t=0:dt:2;                     %Range of time domain

W=2;                           %Width of x-axis
L=2;                           %Length of y-axis
H=2;                           %Height of z-axis

M=30;                          %Number of series expansions for index i
N=30;                          %Number of series expansions for index j
P=30;                          %Number of series expansions for index k

Ti=300;                        %Initial Temperature at t=0 for all x, y, and z
kb=0.003;                      %Heat Conduction Coefficient in W/(m*K)

%%
% Analytical Solution

% Preallocating matrix T(i,j,k,l)
T = zeros(nx,ny,nz,nt);

for l=1:nt
    for i=1:nx
        for j=1:ny
            for k=1:nz
                for m = 1:M
                    for n = 1:N
                        for p = 1:P
                         
        Amnl = ((((2*m-1)*pi)/(W))^2)+((((2*n-1)*pi)/L)^2)+((((2*p-1)*pi)/H)^2);
        mum = ((2*m-1)*pi)/W;
        vun = ((2*n-1)*pi)/L;
        kal = ((2*p-1)*pi)/H;
        
        T(i,j,k,l) = T(i,j,k,l) + ((64*(Ti))/(pi^3))*((sin(mum*x(i))*sin(vun*y(j))*sin(kal*z(k))*exp(-Amnl*kb*t(l)))/((2*m-1)*(2*n-1)*(2*p-1)));
  
                        end
                    end
                end
            end
        end
    end
end

%%
% Postprocessing - Plot Results
Tp=T(4,4,4,:);

surf(x,y,Tp,'EdgeColor','none');
shading interp
title({'3D Steady State Laplace Equation - Analytical Solution'})
xlabel('X-axis (x) \rightarrow')
ylabel('{\leftarrow} Y-axis (y)')
zlabel('Temperature Profile (T) \rightarrow')            
            
