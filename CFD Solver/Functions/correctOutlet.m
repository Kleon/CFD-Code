%Purpose: Calculates the outlet correciton velocity ucorr and corrects
%the outlet velocities to enforce global mass conservation 

%Input: 
%u: staggered velocity u*
%v: staggered velocity v* 

%Output: 
%u: staggered velocity u* with potential outlet correction applied
%v: staggered velocity v* with potential outlet correction applied

%Global:
%h:   mesh spacing h that is the same in both the x and y direction 
%yc:  y-coordinate of cell center incl ghost cells

function [u,v] = correctOutlet(u,v)
global h yc 

%Determine size of v:
[~,N] = size(u);
N = N-2;

%Sum of length of outlets:
Lout = .5;

ucorr = (h/Lout) * ( sum(u(1,2:end-1)) - sum(u(end,2:end-1)) + sum(v(2:end-1,1)) - sum(v(2:end-1,end)) ) ;

%% Left wall [Outlet] u(cell centered):
for j = 2:N+1

        %Left wall (outlet):
        if (yc(j) >= .25) && (yc(j) <= .75)     %(0.25 to 0.75)

        %Apply correction: 
        u(1,j) =  u(1,j) - ucorr;
        end

end

end