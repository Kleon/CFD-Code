%Purpose: applies the fractional step method zero-Nuemann boundary
%condition for Φ and the error e of the multigrid method on a cell centered
%mesh with 1 ghost layer

%Input: 
%phi: cell centered mesh solution variable with 1 ghost layer

%Output: 
%phi: cell centered mesh solution variables with boundary conditions
%applied

function [phi] = bcGS(phi)

%Variables: 
[M, N] = size(phi);          % N = Interior cells in the y direction (j)
                             % M = Interior cells in the x direction (i)

%Subtract out the ghost cell values such that numnber of interior nodes are left: 
N = N-2;
M = M-2;


for j = 1:N+1    %Loop through interior nodes in y direction 
     for i = 1:M+1  %Loop through interior nodes in x direction
    
        phi(1,j) = phi(2,j);        %Left    
        phi(end,j) = phi(end-1,j);  %Right
        phi(i,1) = phi(i,2);        %Bottom
        phi(i,end) = phi(i,end-1);  %Top
        
     end
end
end
