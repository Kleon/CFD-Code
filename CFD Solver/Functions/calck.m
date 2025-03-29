%Purpose: calculates k(t) using the composite 2D midpoint integration rule

%Input: 
%u: 2D array of staggered mesh velocities u 
%v: 2D array of staggered mesh velocities v

%Output: 
% k: kinetic energy in chamber [kx+ky] from u and v velocities 

%Global: 
%h: mesh spacing

function [k] = calck(u,v)
global h 
%% x-direction [kx]%%

%Determine size of u: 
[Mi,Ni] = size(u);
M = Mi-1;    %Interior nodes (X) [node based]
N = Ni-2;    %Interior nodes (Y) [cell centered]

%Preallocate matrix du:
uhatx = zeros(Mi+1,Ni);        %cell centered mesh with 1 ghost layer  
uhaty = zeros(Mi+1,Ni);        %cell centered mesh with 1 ghost layer  

%Calculate kx utilizing 2D midpoint integration rule: 
for j = 2:N+1       %Loop over y
    for ii  = 1:M   %Loop over x

    uhatx(ii+1,j) = ((u(ii+1,j) + u(ii,j)) / (2))^2 ; 

    end
end

%% y-direction [ky]%% 
%Determine size of v:

[Mi,Ni] = size(v);
M = Mi-2;    %Interior nodes (X) [cell centered]
N = Ni-1;    %Interior nodes (Y) [node based]

%Calculate kx utilizing 2D midpoint integration rule: 
for ii  = 2:M+1       %Loop over x
    for j = 1:N       %Loop over y

    uhaty(ii,j+1) =  ((v(ii,j+1) + v(ii,j)) / (2))^2 ; 

    end
end

%% Combine kx and ky to get total kinetic energy in chamber %%

uhat = (1/2)* (uhatx+uhaty)*h^2;

k = sum(sum(uhat));



end