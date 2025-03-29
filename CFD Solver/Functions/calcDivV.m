%Purpose: Calculates the divergence of a staggered velocity field ∇*V in
%the interior of a cell centered mesh with 1 ghost layer. Set the
%divergence in the ghost layers to 0. 

%Input:
%u: staggered velocity u* with outlet correction applied 
%v: staggered velocity v* with outlet correciton applied 

%Output: 
%divV: velocity divergence on a cell centered mesh with 1 ghost layer

%Global: 
%h: mesh spacing 

function [divV] = calcDivV(u,v)
global h 

%Determine size of u: 
[Mi,Ni] = size(u);
M = Mi-1;    %Interior nodes (X) node based
N = Ni-2;    %Interior nodes (Y) cell centered

%Preallocate matrix divV:
divV = zeros(Mi+1,Ni);        %cell centered mesh with 1 ghost layer    


%Calculate du/dx in interior: 
for j = 2:N+1       %Loop over y
    for ii  = 1:M   %Loop over x

    divV(ii+1,j) = (u(ii+1,j) - u(ii,j)) / (h) ;  %Calculate [du/dx]

    end
end

%Determine size of V: 
[Mi,Ni] = size(v);
M = Mi-2;    %Interior nodes (X) cell centered 
N = Ni-1;    %Interior nodes (Y) node based

%Calculate [dv/dy] in interior and add to [du/dx] to calculate the divergence: 
for j = 1:N       %Loop over y
    for ii  = 2:M+1   %Loop over x
        
    dv =  (v(ii,j+1) - v(ii,j)) / (h);  %Calculate [dv/dy]
    divV(ii,j+1) = divV(ii,j+1) + dv ;  % calculate the divergence: div(v) = [du/dx] + [dv/dy]

    end
end


end