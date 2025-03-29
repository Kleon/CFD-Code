%Purpose: Calculates the hyperbolic term Hy of Eq. (3) [highlighted in red] 
% when moved to the right hand side, using the FTCS scheme with 1 ghost layer. Inside the
%function use second order accuracy to calculate the velocities at cell
%centers. 

%Input:
%Y:  Mass fraction Y on a cell centered mesh with 1 ghost layer at a time t^n
%u:  u velocity on a staggered mesh at time t^n
%v:  v velocity on a staggered mesh at time t^n

%Output: 
%HY: hyperbolic terms of HY of Eq. (3) on a cell centered mesh with 1 ghost
%layer

%Global variables: 
%h: mesh spacing

function [HY] = hyperbolic_Y_FTCS_2D(Y,u,v)
global h

%Define M and N: 
[Mi,Ni] = size(Y);
M = Mi-2;           %# Interior Nodes in (x) - Cell centered
N = Ni-2;           %# Interior Nodes in (y) - Node based

%Preallocate variable:
HY = zeros(size(Y));

%Begin loop:
for j = 2:N+1          %Loop through interior (y)
    for ii = 2:M+1       %Loop through interior (x)

    HY(ii,j) = ( 1/2*(u(ii-1,j) + u(ii,j)) )* ( (Y(ii+1,j)-Y(ii-1,j)) / (2*h) );
    
    HY(ii,j) = HY(ii,j) + ...
        (1/2* (v(ii,j-1)+v(ii,j)) ) *...
        ( (Y(ii,j+1)-Y(ii,j-1)) / (2*h) );
    
    end
end

%Account for RHS: 
HY = -HY;

end