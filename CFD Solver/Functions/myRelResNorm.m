function [rr] = myRelResNorm(phi,f,h)
%This function will calculate the infinity norm of a relative residual for
%the two dimensional Poisson equation

%Input: 
%Phi: 2D array of cell centered PDE solution variable including ghost cells
%f: 2D array of cell centered PDE right hand side values including ghost cells
%h: Mesh spacing (equidstant mesh spacing in the x and y direction) 

%Output: 
%rr: infinity norm of relative resiudal

%% Calculate the code %%
%Variables: 
[M, N] = size(phi);          % N = Interior cells in the y direction (j)
                             % M = Interior cells in the x direction (i) 

%Subtract out the ghost cell values such that numnber of interior nodes are left: 
N = N-2;
M = M-2;

%Calculate the residuals:
[r] = myResidual(phi,f,h);
%Calculate the Rinf of the residuals: 
Rinf = max(max(abs(r(2:M+1,2:N+1))));

%Calculate the finf of the solution: 
finf = max(max(abs(f(2:M+1,2:N+1))));

%Calculate the relative residual [rr]:
rr = Rinf/finf;

end