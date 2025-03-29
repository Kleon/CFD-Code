% AEE471 - Final Exam 

%Purpose: This function calculates the time step Δt defined as Δt =
% CFL*Δtmax potentially adjusted to reach a requested output time
% t_out. Here Δt_max is the maximum stable time step of the system of Eqs
% 1-3 when solved using the the FTCS method. 

%Inputs:
% t:                current time t^n
% outputTime:       next requested time for output t_out
%u:                 u velocity on a staggered mesh at t^n
%v:                 v velocity on a staggered mesh at t^n

%Outputs:
% dt:               time step Δt to perform
% outputFlag:       flag, set to 1 if dt was adjusted to reach outputTime, else
% set to 0

%Global Variables: 
% Re: Reynolds number
% Sc: Schmidt number 
%  h: mesh spacing
%CFL: security factor of Δt

function [dt, outputFlag] = calcDtFinal471(t, outputTime,u,v)
    %global Re h CFL 
    global  h CFL 

%% Calculate Δt hyperbolic %%
    %Since calcDtFinal2471 indicated that dthyperbolic is always smaller
    %will only utilize hyperbolic time step for speed:

    dt = CFL * (h/ ( max(max(abs(u))) + max(max(abs(v)))));

    %Calculate hyperbolic dt's (u,v,Y):
    %dt_hyperbolic = CFL * (h/ ( max(max(abs(u))) + max(max(abs(v)))));


 %% Calculate Δt parabolic %%
    % %Assign the Nu for (u,v,Y): 
    % Nu = (1/Re);
    % 
    % %Calculate parabolic dt's (u,v,Y):
    % dt_parabolic = CFL * ((h^2*h^2)/(2*Nu* (h^2+h^2)));

 %% determine Δt to use %%
%Pick the smallest value as Δtmax:
% dt = min(dt_hyperbolic,dt_parabolic);

%% Determine Flags
    if (t < outputTime) && (t + dt >= outputTime)

        dt = outputTime - t;

        %At the time of interest:
        outputFlag = 1;
    else

        %Not at time of interest: 
        outputFlag = 0;

    end

end
