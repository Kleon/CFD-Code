% AEE471 - Exam 7 - Problem 3

%Purpose: This function calculates the hyperbolic terms of the Burgers eqs. (1) & (2)
%[highlighted in red] using second-order central differences on a staggered
%mesh 

%Inputs:
%u:                 u velocity on a staggered mesh at t^n
%v:                 v velocity on a staggered mesh at t^n

%Outputs:
%Hu: hyperbolic terms of Eq. 1 on a staggered mesh 
%Hv: hyperbolic terms of Eq. 2 on a staggered mesh 

%Global Variables: 
%h: mesh spacing 

function [Hu,Hv] = hyperbolic_uv_2D(u,v)
    global h

%% Preallocation & defining meshes %%

   %Preallocate Variables: 
   Hv      = zeros(size(v));
   Hu      = zeros(size(u));

   %Determine size of mesh [u]: 
   [Mu,Nu] = size(u);
   Mu = Mu-2;         %M: Interior points in x [cell centered]
   Nu = Nu-1;         %N: Interior points in y [Node]

   %Determine size of mesh [v]: 
    [Mv,Nv] = size(v);
    Mv = Mv-1;         %M: Interior points in x [cell centered]
    Nv = Nv-2;         %N: Interior points in y [Node]

   %% Calculating Hyperbolic terms of EQ.1 [Hu]: [u] %%

   for j = 2:Nu
       for ii = 2:Mu+1   
       
       %Calculate [duu/dx]
       Hu(ii,j) = ( ((1/2) * (u(ii,j) + u(ii+1,j))) ^2 - ...
            ((1/2) * (u(ii,j) + u(ii-1,j))) ^2 ) / h; 
       
       %Calculate [duv/dy]
       Hu(ii,j) = Hu(ii,j) + ...
            ( ( ( (1/2) * ((u(ii,j)) + u(ii,j+1))) * ((1/2) * ((v(ii,j)) + v(ii+1,j))) ) ...
            - ( ((1/2) * ((u(ii,j)) + u(ii,j-1)))   *  ((1/2) * ((v(ii,j-1)) + v(ii+1,j-1))) ) ) / h; 

       end
   end
   
%Account for RHS: 
Hu = -Hu;

   %% Calculating Hyperbolic terms of EQ.2 [Hv]: [v] %%

   for j = 2:Nv+1
       for ii = 2:Mv

       %Calculate [dvv/dy]
       Hv(ii,j) = ( ((1/2) * (v(ii,j) + v(ii,j+1))) ^2 - ...
            ((1/2) * (v(ii,j) + v(ii,j-1))) ^2 ) / h; 
       
       %Calculate [duv/dy]
       Hv(ii,j) = Hv(ii,j) + ...
            ( ( ( (1/2) * ((u(ii,j)) + u(ii,j+1))) * ((1/2) * ((v(ii,j)) + v(ii+1,j))) ) ...
            - ( ((1/2) * ((u(ii-1,j)) + u(ii-1,j+1)))   *  ((1/2) * ((v(ii-1,j)) + v(ii,j))) ) ) / h; 

       end
   end
%Account for RHS: 
Hv = -Hv;

end



