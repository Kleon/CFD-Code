function [r2h] = myRestrict(rh)
%input:
%rh  = residual on a fine cell centered mesh 

%Output:
%r2h = restricted residual on a coarse cell centered mesh 

%% Function code %%

%Determine interior nodes M (x) and N (y):
[M,N]=size(rh);   

%Subtract out the ghost cell values such that number of interior nodes are left: 
N = N-2;
M = M-2;

%Mh -> M2h (fine to coarse) [4 nodes are go to 1 node]
M = (M/2);
N = (N/2);

% %Preallocate variables: 
r2h = zeros(M+2,N+2);

%Loop through residual matrix and apply restriction (fine to coarse)
   for j = 2:N+1    %Loop through interior nodes in y direction 
        for i = 2:M+1  %Loop through interior nodes in x direction
        
         %Apply restriction: 
         r2h(i,j) = .25*( rh((2*i)-2,(2*j)-2) + rh((2*i)-1,(2*j)-2) +rh((2*i)-2,(2*j)-1) + rh((2*i)-1,(2*j)-1) );

        end
   end

end