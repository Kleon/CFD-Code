function [v] = parabolic_FTCS_2D_v(v,Qu,dt)
%Purpose: Performs single time step of the FTCS method for Eq 2 on a
%staggered mesh for Exam5

%Input:
%v:  2D array of staggered mesh values of v^n
%Qu: 2D array of staggered mesh values of Qu at t^n
%dt: time step to perform deltaT

%Output: 
%v: 2D array of staggered mesh values v^n+1 with boundary conditions
%applied

%Global variables: 
%Re: value of the reynolds number
%t:  current time t^n
%h:  mesh spacing

%% Code %% 
global Re t h

%Define M and N: 
[Mi,Ni] = size(v);
M = Mi-2;           %# Interior Nodes in (x) - Cell centered
N = Ni-1;           %# Interior Nodes in (y) - Node based

%Preallocate variable:
unew = zeros(size(v));

%Begin loop:
for j = 2:N          %Loop through interior (y)
    for ii = 2:M+1       %Loop through interior (x)
    unew(ii,j) = (dt/(Re*h^2)) * ( (v(ii+1,j)-2*v(ii,j)+v(ii-1,j)) + (v(ii,j+1)-2*v(ii,j)+v(ii,j-1)) ) + v(ii,j) + dt*Qu(ii,j);    %T: Temp at next time step
    end
end

%Set boundary conditions 
v = bc_v(unew,t+dt);


end