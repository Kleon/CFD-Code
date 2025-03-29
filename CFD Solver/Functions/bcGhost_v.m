%Purpose: applies only the ghost cell boundary conditions for the v
%velocity on a staggered mesh with at least second order accuracy in space

%Input: 
%u: 2D array of staggered mesh values of v
%t: time at which to apply the boundary conditions


%Output: 
%v: 2D array of staggered mesh values of v with ghost cell boundary
%conditions applied 

%Global: 
%yf: y-coordinates of cell faces


function [v] = bcGhost_v(v,t)
global yf

%Ghost cells for u are located at top and bottom. Left are right can be left unchanged: 

%Define givens from problem:
alpha1 = pi;

U1 = 3/4;

f = @(z) 6*z*(1-z); 

%Define the no slip condition: 

%Dirichlet BC: [g = 0] no slip, cell centered 
v(end,:) = -v(end-1,:);   %Right boundary condition
v(1,:)   = -v(2,:);       %Left boundary condition

%Determine the size of v: 
[~, N] = size(v);
N = N-1;                %Interior cells in x

%% Left wall: 
for j = 1:N+1
        %left wall (outlet + no slip) & [cell centered]:  
        if (yf(j) >= .25) && (yf(j) <= .75)         %(.25 to .5)

        %Nuemann BC
        v(1,j) = v(2,j);
        end
end

%% Right wall:
for j = 1:N+1
        %Right wall (Inlet1 + no slip) & [cell centered]: 
        if (yf(j) >=.25) && (yf(j) <=1.25)     %(.25 to 1.25)
        z1 = (yf(j) - 0.25)/1;
        vin1 = (U1*sin(alpha1)*f(z1));         %Inlet 1

        %Dirichlet BC
        v(end,j) = 2*vin1 -v(end-1,j);
        end
end


end