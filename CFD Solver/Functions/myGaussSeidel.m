function [phi] = myGaussSeidel(phi,f,h,niter)
%Input: 
%phi: 2D array of cell centered PDE solution variable including ghost cells
%f: array of cell centered PDE right hand side values including ghost cells
%h: mesh spacing (equidistant mesh spacing in x and y direction)
%niter: integer number of Gauss Seidel iterations to be performed

%Output:
%Phi: 2D array of cell centered PDE solution variable including ghost
%cells, containing the solution after niter iterations of Gauss Seidel

%% Code %%
%Working with cell centered mesh 
[M, N] = size(phi);          % N = Cells in the x direction (i)
                             % M = Cells in the y direction (j) 

%Subtract out the ghost cell values such that number of interior nodes are left: 
N = N-2;
M = M-2;

%Gauss-Seidel Loop: 
for k=1:niter       %Loop through requested iterations

   for j = 2:N+1    %Loop through interior nodes in y direction 
        for i = 2:M+1  %Loop through interior nodes in x direction
        
         %Apply Gauss Seidel :
         phi(i,j)= .25*( phi(i-1,j) + phi(i+1,j) + phi(i,j-1) + phi(i,j+1)) - (.25*h^2*f(i,j)); 

        end
   end

%Provides boundary conditions:
phi = bcGS(phi);
end






end