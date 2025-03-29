function [Y] = parabolic_FTCS_2D_Y(Y,Qu,dt)

%Purpose: Performs single time step of the FTCS method for Eq 3 on a
%cell centered mesh for Exam5

%Input:
%Y:  2D array of staggered mesh values of Y^n
%Qu: 2D array of staggered mesh values of QY at t^n
%dt: time step to perform deltaT

%Output: 
%Y: 2D array of staggered mesh values v^n+1 with boundary conditions
%applied

%Global variables: 
%Re: value of the reynolds number
%Sc: value of tthe schmidt number
%t:  current time t^n
%h:  mesh spacing

%% Code %% 
global Sc Re t h

%Define M and N: 
[Mi,Ni] = size(Y);
M = Mi-2;           %# Interior Nodes in (x) - Cell centered
N = Ni-1;           %# Interior Nodes in (y) - Node based

%Preallocate variable:
unew = zeros(size(Y));

%Begin loop:
for j = 2:N          %Loop through interior (y)
    for ii = 2:M+1       %Loop through interior (x)
    unew(ii,j) = (dt/(Re*Sc*h^2)) * ( (Y(ii+1,j)-2*Y(ii,j)+Y(ii-1,j)) + (Y(ii,j+1)-2*Y(ii,j)+Y(ii,j-1)) ) + Y(ii,j) + dt*Qu(ii,j);    %T: Temp at next time step
    end
end

%Set boundary conditions 
Y = bc_Y(unew,t+dt);

end