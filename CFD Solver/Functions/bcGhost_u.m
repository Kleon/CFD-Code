%Purpose: applies only the ghost cell boundary conditions for the u
%velocity on a staggered mesh with at least second order accuracy in space

%Input: 
%u: 2D array of staggered mesh values of u 
%t: time at which to apply the boundary conditions


%Output: 
%u: 2D array of staggered mesh values of u with ghost cell boundary
%conditions applied 

%Global: 
%xf: x-coordinates of cell faces


function [u] = bcGhost_u(u,t)
global xf

%Ghost cells for u are located at top and bottom. Left are right can be left unchanged: 

%Define givens from problem:
alpha2 = pi/3;
alpha3 = -pi/6;

U2 = 2;
U3 = 3;

f = @(z) 6*z*(1-z); 

%Define the no slip condition: 

%Dirichlet BC: [g = 0] no slip, cell centered
u(:,end) = -u(:,end-1);         %Top boundary condition
u(:,1)   = -u(:,2);             %Bottom boundary condition

%Determine the size of u: 
[M,~] = size(u);
M = M-1;                %Interior cells in x

%% Top wall: 
for i = 1:M+1
        %Top wall (Inlet3 + no slip) & [cell centered]:
        if (xf(i) >= 2) && (xf(i) <= 2.5)     %(2 to 2.5)
        z3 = (xf(i)-2)/.5;
        vin3 = (U3*cos(alpha3)*f(z3));          %Inlet 3

        %Dirichlet BC [g=vin3]
        u(i,end) = 2*vin3 - u(i,end-1);

        end
end

%% Bottom wall:
for i = 1:M+1
        %Bottom wall (Inlet2 + no slip) & [cell centered]:
        if (xf(i) >= 1.5) && (xf(i) <= 2)     %(1.5 to 2)
        z2 = (xf(i)-1.5)/.5;
        vin2 = (U2*cos(alpha2)*f(z2));          %Inlet 2

        %Dirichlet BC [g = vin2]
        u(i,1) = 2*vin2 - u(i,2);
        end
end



end