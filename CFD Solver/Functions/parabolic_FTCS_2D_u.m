function [u] = parabolic_FTCS_2D_u(u,Qu,dt)

%Input:
%u:  2D array of staggered mesh values of u^n
%Qu: 2D array of staggered mesh values of Qu at t^n
%dt: time step to perform deltaT

%Output: 
%u: 2D array of staggered mesh values 

%Global variables: 
%Re: value of the reynolds number
%t:  current time t^n
%h:  mesh spacing

%% Code %% 
global Re t h

%Define M and N: 
[Mi,Ni] = size(u);
M = Mi-1;           %# Interior Nodes in (x) - Node based
N = Ni-2;           %# Interior Nodes in (y) - Cell centered

%Preallocate variable:
unew = zeros(size(u));

%Begin loop:
for j = 2:N+1          %Loop through interior (y)
    for ii = 2:M       %Loop through interior (x)
    unew(ii,j) = (dt/(Re*h^2)) * ( (u(ii+1,j)-2*u(ii,j)+u(ii-1,j)) + (u(ii,j+1)-2*u(ii,j)+u(ii,j-1)) ) + u(ii,j) + dt*Qu(ii,j);    %T: Temp at next time step
    end
end

%Set boundary conditions 
u = bc_u(unew,t+dt);


end