function [u] = bc_u(u,t)
%Purpose: This funciton will apply boundary conditions for the u velocity
%on a staggered mesh with at least second order accuracy in space

%Input: 
%u: 2D array of staggered mesh values of u 
%t: time at which to apply the boundary conditions

%Output: 
%u: 2D array of staggered mesh values of u with BCs applied

%Global variables: 
%xf: x-coordinate of cell faces
%yc: y-coordinate of cell center incl ghost cells

%% Code %%
global xf yc 
%Define givens from problem:
alpha1 = pi;
alpha2 = pi/3;
alpha3 = -pi/6;
U1 = 3/4;
U2 = 2;
U3 = 3;
f = @(z) 6*z*(1-z); 

%Define M and N: 
[Mi,Ni] = size(u);
M = Mi-1;           %# Interior Nodes in (x) - Node based
N = Ni-2;           %# Interior Nodes in (y) - Cell centered

%Define the no slip condition on chamber walls: 
u(end,:) = 0;               %Right boundary condition
u(1,:)   = 0;               %Left boundary condition

u(:,N+2) = -u(:,N+1);       %Top boundary condition
u(:,1)   = -u(:,2);         %Bottom boundary condition

%To prevent intereference of intervals seperate for loops are used:

%% Left wall:
for j = 1:N+2
        %left wall (outlet + no slip) & [Node based]:  
        if (yc(j) >= .25) && (yc(j) <= .75)         %(.25 to .5)
        %Nuemann BC
        u(1,j) = 1/3*(4* u(2,j) - u(3,j));
        end
end

%% Right wall:
for j = 1:N+2
        %Right wall (Inlet1 + no slip) & [Node based]: 
        if (yc(j) >=.25) && (yc(j) <=1.25)     %(.25 to 1.25)
        z1 = (yc(j)-.25)/1;
        %Dirichlet BC
        u(end,j) = U1*cos(alpha1)*f(z1); 
        end
end

%% Top wall: 
for i = 1:M+1
        %Top wall (Inlet3 + no slip) & [cell centered]:
        if (xf(i) >= 2) && (xf(i) <= 2.5)     %(2 to 2.5)
        z3 = (xf(i)-2)/.5;
        %Dirichlet BC
        u(i,end) = 2*(U3*cos(alpha3)*f(z3)) - u(i,end-1);
        end
end

%% Bottom wall:
for i = 1:M+1
        %Bottom wall (Inlet2 + no slip) & [cell centered]:
        if (xf(i) >= 1.5) && (xf(i) <= 2)     %(1.5 to 2)
        z2 = (xf(i)-1.5)/.5;
        %Dirichlet BC
        u(i,1) = 2*(U2*cos(alpha2)*f(z2)) - u(i,2);
        end
end


end