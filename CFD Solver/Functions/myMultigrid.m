function [phi] = myMultigrid(phi,f,h)
%Purpose: The function shall perform a single V-cycle using a single Gauss-Seidel 
% iteration on each mesh level with an additional 3 iterations on the coarsest mesh 
% (nc=3 in the slides), and output phi, the 2D array of the cell centered 
% solution variable including updated ghost cells.

%Input:
%Phi = 2D array of cell centered solution variable including ghost cells
%f   = 2D array of cell centered PDE right hand side values including ghost
%cells 

%h   = equidistant mesh spacing in the x and y direction 

%Output:
%phi = 2D array of cell centered solution variable including updated ghost
%cells 


%% Code %%

%Define interior nodes for M & N (x & y)
M = size(phi,1)-2;
N = size(phi,2)-2;

%Run initial GaussSeidel:
n = 1;
phi = myGaussSeidel(phi,f,h,n);

if mod(M,2) == 0 && mod(N,2) == 0   %Runs if M and N are divisable by 2 
    rh  = myResidual(phi,f,h);      %Calculate residual
    r2h = myRestrict(rh);           %restrict residual to coarses grid
    e2h = zeros(M/2+2,N/2+2);       %2D cell centered
    e2h = myMultigrid(e2h,r2h,2*h); %solve course mesh 
    eh = myProlong(e2h);            %Prolong error to fine grid
    phi = phi + eh;                 %apply correction to grid solution
    phi = bcGS(phi);                %apply correction to grid solution
    phi = myGaussSeidel(phi,f,h,n);%On fine grid perform n iterations 
else

    %Perform 3 iterations of Gauss Seidel on coarses mesh
    nc=3;
    phi = myGaussSeidel(phi,f,h,nc);
end


end