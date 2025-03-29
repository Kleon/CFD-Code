%% AEE471 - Exam - Problem 4

%Purpose: Applies the boundary conditions for the v velocity on a staggered
%mesh with at least second order accuracy in space

%Inputs:
% v: 2D array of staggered mesh values of v
% t: time at which to apply the boundary conditions

%Outputs:
% v: 2D array of staggered mesh values of v with boundary conditions
% applied

%Global Variables: 
% xf: x-coordinates of cell centers incl ghost cells
% yf: y-coordinates of cell faces

%% Code %%
function [v] = bc_v(v, t)
global xc yf 

%Define givens from problem:
alpha1 = pi;
alpha2 = pi/3;
alpha3 = -pi/6;

U1 = 3/4;
U2 = 2;
U3 = 3;

f = @(z) 6*z*(1-z); 

%Define M and N: 
[Mi,Ni] = size(v);
M = Mi-2;           %# Interior Nodes in (x) - Node based
N = Ni-1;           %# Interior Nodes in (y) - Cell centered

%Define the no slip condition on chamber walls: 
v(end,:) = -v(end-1,:);   %Right boundary condition
v(1,:)   = -v(2,:);   %Left boundary condition
v(:,end) = 0;         %Top boundary condition
v(:,1)   = 0;         %Bottom boundary condition

%To prevent intereference of intervals seperate for loops are used:

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

%% Top wall: 
for i = 1:M+2
        %Top wall (Inlet3 + no slip) & [node based]:
        if (xc(i) >= 2) && (xc(i) <= 2.5)     %(2 to 2.5)
        z3 = (xc(i)-2)/.5;
        vin3 = (U3*sin(alpha3)*f(z3));          %Inlet 3

        %Dirichlet BC
        v(i,end) = vin3;
        end
end

%% Bottom wall:
for i = 1:M+2
        %Bottom wall (Inlet2 + no slip) & [node based]:
        if (xc(i) >= 1.5) && (xc(i) <= 2)     %(1.5 to 2)
        z2 = (xc(i)-1.5)/.5;
        vin2 = (U2*sin(alpha2)*f(z2));          %Inlet 2

        %Dirichlet BC
        v(i,1) = vin2;
        end
end


end