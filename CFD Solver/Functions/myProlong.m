function [eh] = myProlong(e2h)
%Purpose: This function will take as input the error (e2h) on a coarse cell
%centered mesh and return the prolonged error (eh) on a fine cell
%centeredmesh. 

%Input:
%e2h = error on a coarse cell centered mesh

%Output:
%eh  = prolonged error on a fine cell centered mesh

%% Code %%

%Determine interior nodes M (x) and N (y):
[M,N]=size(e2h);   

%Subtract out the ghost cell values such that number of interior nodes are left: 
N = N-2;
M = M-2;

%M2h -> Mh (coarse to fine) [1 nodes go to 2 node]
Mh = (M*2);
Nh = (N*2);

%Preallocate variables: 
eh = zeros(Mh+2,Nh+2);       %Include the ghost cells

%Loop through residual matrix and apply prolongation (coarse to fine)
   for j = 2:N+1    %Loop through interior nodes in y direction 
        for i = 2:M+1  %Loop through interior nodes in x direction
        
         %Apply prolongation: [1 node to 4]
         eh((2*i)-2,(2*j)-2) = e2h(i,j);
         eh((2*i)-1,(2*j)-2) = e2h(i,j);
         eh((2*i)-2,(2*j)-1) = e2h(i,j);
         eh((2*i)-1,(2*j)-1) = e2h(i,j);
        
        end
   end

%Provides boundary conditions:
eh = bcGS(eh);

end