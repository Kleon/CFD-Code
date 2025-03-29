function [phi,Linf,iter] = myPoisson(phi,f,h,nIterMax,epsilon)
%Purpose: 
% This function shall perform at least 1 V-cycle iteration and
%continue to perform V-cycle iterations until either the maximum number of
%iterations set by nIterMax has been reached or infinity norm of the
%relative residual is smaller than epsilon

%Input:
%Phi = 2D array of cell centered solution variable including ghost cells
%f   = 2D array of cell centered PDE right hand side values including ghost
%cells 
%h   = mesh spacing (the equidistant mesh spacing in the x and y direction)
%epsilon = convergence threshold for infinity norm of the relative residual

%output:
%Phi = 2D array of cell centered solution variable including ghost cells 
%Linf = Infinity norm of the relative residual of the solution
%epsilon = convergence threshold for infinity norm of the relative resiudal

%% V Cycle Code %%

%Define the logical conditions to compute the multigrid method: 
for iter = 1:nIterMax     %Compute the number of iterations defined by nIterMax

        [phi] = myMultigrid(phi,f,h); %Compute the V method
        [Linf] = myRelResNorm(phi,f,h);  %Compute the infinity norm of the residual

    if Linf < epsilon  || iter == nIterMax             %Check if infinity norm is smaller than epsilon
        break           %End if infinity norm norm of relative residual is smaller than epsilon
    end
end