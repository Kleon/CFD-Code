%Purpose: Projects the staggered velocity v into the subspace of divergence
%free velocity fields. After projection apply an update to the ghost cell
%velocities t^n+1. Do not apply an update to velocities directly on the
%boundaries. 

%Input: 
%u:     staggered velocity u* 
%v:     staggered velocity v* 
%phi:   cell centered lagrange multiplier phi with 1 ghost cell layer 
%dt:    time step deltat

%Output: 
%u:     divergence free staggered velocity u^n+1
%v:     divergence free staggered velocity v^n+1

%Global variables: 
%t:     current time t^n
%h:     mesh spacing h 


function [u,v] = projectV(u,v,phi,dt)
global t h 
%% u [left to right] %%
%Determine the size of staggered mesh: 
[M,N] = size(u);
M = M-1;            %Node based in (x) 
N = N-2;            %Cell centered in (y)

%Calculate u in interior: 
for j = 2:N+1       %Loop over y
    for ii  = 2:M   %Loop over x

    u(ii,j) =  u(ii,j) - dt * ( (phi(ii+1,j) - phi(ii,j)) / (h) );  %project into solenoidal subspace

    end
end

%Update the ghost cell values: 
[u] = bcGhost_u(u,t+dt);

%% v [bottom to top] %%
%Determine the size of staggered mesh: 
[M,N] = size(v);
M = M-2;            %Cell centered in (x) 
N = N-1;            %Node based in (y)

%Calculate v in interior:
for ii  = 2:M+1   %Loop over x
    for j = 2:N   %Loop over y

    v(ii,j) =  v(ii,j) - dt * ( (phi(ii,j+1) - phi(ii,j)) / (h) );  %project into solenoidal subspace

    end
end

%Update the ghost cell values: 
[v] = bcGhost_v(v,t+dt);

end