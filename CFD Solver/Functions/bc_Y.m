function [Y] = bc_Y(Y,t)
%Purpose: Applies the boundary conditions for the mass fraction Y on a cell
%centered mesh with 1 ghost cell layer using second order accuracy in space

%Input:
%Y: 2D array of cell centered values of Y
%t: time at which to apply the boundary conditions 

%Output: 
%Y: 2D array of cell centered values of Y with boundary conditions applied

%Global variables: 
%xc: x-coordinates of cell centers incl. ghost cells
%yx: y-coordinates of cell centers incl. ghost cells


%% Code %%
global xc yc

%Define M and N: 
[Mi,Ni] = size(Y);
M = Mi-2;           %# Interior Nodes in (x) - cell centered
N = Ni-2;           %# Interior Nodes in (y) - Cell centered

%inlet conditions:
Y1 = 1;
Y2 = 0;
Y3 = 0; 

%Fluid in the chamber is initially at rest Y(x,y,t=0) = 0:
Y(end,:) = Y(end-1,:);       %Right boundary condition
Y(1,:)   = Y(2,:);           %Left boundary condition

Y(:,end) = Y(:,end-1);       %Top boundary condition
Y(:,1)   = Y(:,2);           %Bottom boundary condition

%% Left wall:
for j = 1:N+2
        %left wall (outlet + no slip) & [cell centered]:  
        if (yc(j) >= .25) && (yc(j) <= .75)         %(.25 to .5)

        %Nuemann BC
        Y(1,j) = Y(2,j);
        end
end

%% Right wall:
for j = 1:N+2
        %Right wall (Inlet1 + no slip) & [cell centered]: 
        if (yc(j) >=.25) && (yc(j) <=1.25)     %(.25 to 1.25)

        %Dirichlet BC
        Y(end,j) = 2*Y1 - Y(end-1,j);
        end
end

%% Top wall: 
for i = 1:M+2
        %Top wall (Inlet3 + no slip) & [node based]:
        if (xc(i) >= 2) && (xc(i) <= 2.5)     %(2 to 2.5)

        %Dirichlet BC
        Y(i,end) = -Y(i,end-1);
        end
end

%% Bottom wall:
for i = 1:M+2
        %Bottom wall (Inlet2 + no slip) & [node based]:
        if (xc(i) >= 1.5) && (xc(i) <= 2)     %(1.5 to 2)

        %Dirichlet BC
        Y(i,1) = -Y(i,2);
        end
end