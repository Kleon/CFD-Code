%Purpose: calculates eq(11) using the composite 2D midpoint integration
%rule

%Input: 
%Y: 2D array of cell centered values of Y incl ghost cells

%Output: 
%S: value of S

%Global: 
%Lx: size of the domain in the x-direction Lx
%Ly: size of the domain in the x-direciton Ly 
%h : mesh spacing

function [S] = calcS(Y)
global Lx Ly h

%Preallocate S matrix: 
S = zeros(size(Y));

%Cell centered in x & y:
[M,N] = size(Y);
M = M-2;
N = N-2;

%Define loop to solve for 
for ii = 2:M+1
    for j = 2:N+1

    S(ii,j) = (Y(ii,j)* (1 - Y(ii,j))) *h^2;

    end
end

%Perform integration
S = sum(sum(S));
S = (1/(Lx*Ly)).*S;

end