function [r] = myResidual(phi,f,h)
%% Description of the problem: 

%This function will calculate the residual r for the two-dimensional
%poisson equation 

%Inputs:
%Phi:   2D array of cell centered PDE SOLUTIONS variable including ghost cells
%f:     2D array of cell centered PDE RIGHT HAND SIDE values including ghost cells 
%h:     MESH SPACING (the equidistant mesh spacing in the x and y- direction) 

%Output:
%r:     2D array of RESIDUAL calculated on cell centered mesh including ghost cells set to zero

%% Code %% 
%Variables: 
[M, N] = size(phi);          % N = Interior cells in the y direction (j)
                             % M = Interior cells in the x direction (i) 

%Subtract out the ghost cell values such that numnber of interior nodes are left: 
N = N-2;
M = M-2;

%Preallocate variable for loop: 
r = zeros(size(phi));

%Begin the for loop to solve the two-dimensional poisson equation: 

for j = 2:N+1    %Loop through interior nodes in y direction 
    for i = 2:M+1  %Loop through interior nodes in x direction
    
    %Calculate the actual function value: 
    %For interior cells: 
    %Residual (r) = actual (f) - approx (finite differences) 
    r(i,j) = f(i,j) - ( phi(i+1,j) - 4*phi(i,j) + phi(i-1,j) + phi(i,j+1) + phi(i,j-1)) /(h^2);

    end
end

