%Final Project AEE471

%% Problem Statement %%
%Model the flow in the 2D-rectangular mixing chambersize Lx=3 and Ly=4
%using the non-dimensional continuity equation, 2D Navier-Stokesequations
%and a convective/diffusive transport equation with  Sc=1 & Re = 100. The
%chamber contains 2 chemical species absorbing rotating disks that are
%modeled using the Immersed Boundary method covered in class (provided by
%calcSourceIBFinalfunction). The fluid is initially at rest with zero mass
%fraction of chemical species. 

%Determine the flow of: u,v, & Y at desired OutputTimes & performance of
%chamber K & S

clear; clc; 

global xf yc yf xc Re h CFL Sc t Lx Ly

%% Define knowns %% 

%Mesh: 
Lx  = 3;        %Lx: Domain length in x
Ly  = 2;        %Ly: Domain length in y
xs  = 0;        %xs: Starting point of domain in x
ys  = 0;        %ys: Starting point of domain in y 
M   = 96;       %M:Interior nodes in x direction [Node based]
N   = 64;       %N:Interior nodes in y direction [Cell centered]
h  = Ly/N;      %dy: mesh spacing in y direction

%Cell face, cell centered coordinates: 
xf = linspace(xs,xs+Lx,M+1)';            % cell face coordinates
yc = linspace(ys-h/2,ys+Ly+h/2,N+2)';    % cell centered coordinates

yf = linspace(ys,ys+Ly,N+1)';            % cell face coordinates
xc = linspace(xs-h/2,xs+Lx+h/2,M+2)';    % cell centered coordinates

%Givens: 
CFL = .5;                                %CFL: Security factor
Re  = 100;                               %Re: Reynolds
Sc  = 1;                                 %Sc: Schmidt number

%Define times of interest: 
t = 0;                                   %t:  Initial time
tf = 10;                                 %tf: Final time
OutputTimes  = [4,4.25,4.5,4.75,5,5.25,5.5,5.75,6,10];  %Times of interest

%Preallocate variables:
kplot = [];                 %Initialize empty array for storage
Splot = [];                 %Initialize empty array for storage
time  = [];                 %Initialize empty array for storage
dts   = [];                 %Initialize empty array for storage
%timesteptypeStore = [];          %Initialize empty array for storage of which dt is being used to optimize dt function

%% Fractional step initial condition %%

%1) Set initial conditions inside of the computational domain for u,v,Y:

    %Fluid is initially at rest: 
    u = zeros(M+1,N+2);     % cell centered in (y) & cell face in (x)
    v = zeros(M+2,N+1);     % cell centered in (x) & cell face in (y)
    Y = zeros(M+2,N+2);     % cell centered in (x) & (y)
    phi = ones(M+2,N+2);    %cell centered in (x) & (y)

%2) Apply boundary conditions to boundaries and ghost cells using initial time:

    %Ghost cells: 
    [u] = bcGhost_u(u,t);       %Staggered mesh 
    [v] = bcGhost_v(v,t);       %Staggered mesh 

    %Boundaries: 
    [u] = bc_u(u,t);            %Staggered mesh 
    [v] = bc_v(v,t);            %Staggered mesh
    [Y] = bc_Y(Y,t);            %Cell centered mesh 

%3) Correct outlet velocities to ensure volume conservation: 
    [u,v] = correctOutlet(u,v);

%4) Calculate right hand side of Poisson equation using Δt=1:
    dt = 1;                     %Time step
    [divV] = calcDivV(u,v);     %Divergence 
    f = (1/dt)*(divV);          %RHS of equation 

%5) Solve poisson equation for Lagrange multiplier using Neumann boundary conditions:
    nIterMax   = 5;            %Iterations to complete
    epsilon = 1*10^-2;         %Convergence threshold 
    [phi,Linf,iter] = myPoisson(phi,f,h,nIterMax,epsilon);  %Solve
    [phi] = bcGS(phi);         %Apply boundary Conditions

%6) Project velociites using Lagrange multiplier using Δt=1:
%7)Apply boundary conditions to ghost cell velocities only using Initial time

    [u,v] = projectV(u,v,phi,dt);  


%% 8) Start time loop [next portion]

    %k and S for t=0: 
    kplot(1) = calck(u,v);
    Splot(1) = calcS(Y);
    time(1)  = 0;
    dts(1)   = 0;

    %Initialize loop: 
    i=1;            %i:Index for output time
    j=1;            %j:Index for plots
    jj =1;          %jj:  Index for k & S storage 

    OutputTime = OutputTimes(i);    %Intialize Output time to start with

    %Define the figures:
    examFig1 = figure(1);
    examFig2 = figure(2);
    examFig3 = figure(3);
    examFig4 = figure(4);
    examFig5 = figure(5);

    % %Define the location & size of the figures: (in pixels)
    % set(examFig1, 'Position', [50, 500, 500, 500]);   % [left bottom width height]
    % set(examFig2, 'Position', [550, 500, 500, 500]);  % [left bottom width height]
    % set(examFig3, 'Position', [1050, 500, 500, 500]); % [left bottom width height]
    % set(examFig4, 'Position', [1550, 500, 500, 500]); % [left bottom width height]
    % set(examFig5, 'Position', [2050, 500, 500, 500]); % [left bottom width height]


%% Fractional step %% 

%Solve mesh at different discrete times until final time is reached:
while (t<tf)

    %1) Calculate stable time step
        [dt, outputFlag] = calcDtFinal471(t, OutputTime,u,v);

    %Store the timestep type: [calcDtFinal2471]
    %timesteptypeStore(jj) = timesteptype;    %Store the timestep type used

    %2) Solve stokes momentum equations without pressure terms on staggered,
    %and scalar equations (if requested) on a cell centered mesh in interior:
        %Without P. P is used to project u,v back into the subspace of
        %solenoidal functions

    %[Source]
    [Qu,Qv,QY] = calcSourceIBFinal(u,v,Y,t,dt);

    %Hyperbolic: [Convective term]
    [Hu,Hv] = hyperbolic_uv_2D(u,v);
    [HY] = hyperbolic_Y_FTCS_2D(Y,u,v);

    %Parabolic:  [viscous terms]
    [u] = parabolic_FTCS_2D_u(u,Qu+Hu,dt);        %Solve mesh for [u] and apply BC on new time step
    [v] = parabolic_FTCS_2D_v(v,Qv+Hv,dt);        %Solve mesh for [v] and apply BC on new time step
    [Y] = parabolic_FTCS_2D_Y(Y,QY+HY,dt);        %Solve mesh for [Y] and apply BC on new time step  
   
    %3) Apply boundary conditions to boundaries and ghost cells using tn+1
    %Ghost cells: 
    [u] = bcGhost_u(u,t);       %Staggered mesh 
    [v] = bcGhost_v(v,t);       %Staggered mesh 

    %4) Correct outlet velocities to ensure volume conservation
    [u,v] = correctOutlet(u,v);        

    %5) Calculate right hand side of poisson equation [Solutions] (check that sum of rhs =
    %machine zero)
    [divV] = calcDivV(u,v);     %Divergence 
    f = (1/dt)*(divV);          %RHS of equation 

    %6) Solve poisson equation for Lagrange multiplier using Nuemann boundary
    %conditions:
    nIterMax   = 5;            %Iterations to complete
    epsilon = 1*10^-2;      %Convergence threshold 
    [phi,Linf,iter] = myPoisson(phi,f,h,nIterMax,epsilon);

    %7) Project velocities [back to divergence free subspace] using Lagrange multiplier:
    %8) Apply boundary conditions to ghost cell velocities only using tn+1: 

    [u,v] = projectV(u,v,phi,dt);

    %9) Update time and go back to step 1 if output time of interest not
    %met
     t  = t + dt;                              %Update current time
    
    %Calculate k(t) & S(t) for each time step:
    jj = jj+1;                          %Index into next storage place for k(t) & S(t)
    kplot(jj) = calck(u,v);             %Calculate kinetic energy in chamber at time of interest
    Splot(jj) = calcS(Y);               %Calculate the scalar mixing parameter
    time(jj) = t;                       %Store the t value (tnew = t+dt)   
    dts(jj) = dt;                       %Store the dt value

    if outputFlag == 1

        %Verification code is still running: 
        disp(t)

        %Create each subplot: [u,v,Y]
        if t >= 4 && t <= 6

        %u:
            figure(examFig1);
            subplot(3,3,j)
            pcolor(xf, yc, u');      %Create the pseudocolor plot x an y being mesh, z being the u at specific point
            shading interp;  % Interpolated shading
            colormap jet;    % Colormap
            colorbar;        % Colorbar
            % Set title for the subplot
            title(['u at t = ', num2str(OutputTimes(i))]);  %Create dynamic title using outputtime names
            %Labels & limits:
            xlabel('x')         %Labels the x axis
            ylabel('y')         %Labels the y axis
            xlim([0 3])         %Modifies the x limit
            ylim([0 2])         %Modifies the y limit
            clim([-2.5 2.5])    %Modifies the colorbar limits
            %pbaspect([3 2 1])   %Modifies the figure shape to be rectangular
            daspect([2.75 3 1])


        %v:
            figure(examFig2);
            subplot(3,3,j)
            pcolor(xc, yf, v');      %Create the pseudocolor plot x an y being mesh, z being the u at specific point
            shading interp;  % Interpolated shading
            colormap jet;    % Colormap
            colorbar;        % Colorbar
            % Set title for the subplot
            title(['v at t = ', num2str(OutputTimes(i))]);  %Create dynamic title using outputtime names
            %Labels & limits:
            xlabel('x')         %Labels the x axis
            ylabel('y')         %Labels the y axis
            xlim([0 3])         %Modifies the x limit
            ylim([0 2])         %Modifies the y limit
            clim([-2.5 2.5])    %Modifies the colorbar limits
            pbaspect([3 2 1])   %Modifies the figure shape to be rectangular

        %Y:
            figure(examFig3);
            subplot(3,3,j)
            pcolor(xc, yc, Y');      %Create the pseudocolor plot x an y being mesh, z being the u at specific point
            shading interp;  % Interpolated shading
            colormap hot;    % Colormap
            colorbar;        % Colorbar
            % Set title for the subplot
            title(['Y at t = ', num2str(OutputTimes(i))]);  %Create dynamic title using outputtime names
            %Labels & limits:
            xlabel('x')         %Labels the x axis
            ylabel('y')         %Labels the y axis
            xlim([0 3])         %Modifies the x limit
            ylim([0 2])         %Modifies the y limit
            clim([0 1])         %Modifies the colorbar limits
            pbaspect([3 2 1])   %Modifies the figure shape to be rectangular

        %Subplot indexing: 
        j = j+1;                                %index into next spot for plot

        end

        %Output time indexing: 
        i = i+1;                                %index into next output time 

        if i > length(OutputTimes)              %Check if i is greater than time index
            
            %Plots for k & S:
            %k
                figure(examFig4);
                plot(time,kplot);      %Plot output times and kinetic energy at these times
                %Labels & limits:
                title('k vs t')
                xlabel('t')         %Labels the x axis
                ylabel('k')         %Labels the y axis
                xlim([0 10])        %Modifies the x limit

            %S
                figure(examFig5);
                plot(time,Splot);      %Plot output times and scalar mixing parameter at these times
                %Labels & limits:
                title('S vs t')
                xlabel('t')         %Labels the x axis
                ylabel('S')         %Labels the y axis
                xlim([0 10])        %Modifies the x limit
                ylim([0 .25])       %Modifies the y limit

            break                               %Break loop 

            else

        OutputTime = OutputTimes(i);            %get the next output time of interest

        
        end
    end
end

%% Save the figures to .jpg for report %%

%Saving: 
saveas(examFig1, 'uvel.png');
saveas(examFig2, 'vvel.png');
saveas(examFig3, 'Y.png');
saveas(examFig4, 'K_Bar.png');
saveas(examFig5, 'S_Bar.png');



%% Solution Verification Study %%

%Calculate Kbar: 
[K] = myTrapezoidal(time,kplot);
Kbar = (1/10)*K;

%Calculate Sbar: 
[S] = myTrapezoidal(time,Splot);
Rbar = (1/10)*S;

%Create table for excel sheet verification study: 
TK = table(N,Kbar, 'VariableNames', {'N','Kbar'});
TR = table(N,Rbar, 'VariableNames', {'N','Rbar'});


