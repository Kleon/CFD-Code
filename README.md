# CFD-Code

Code to solve Navier-Stokes in 2D in the incompressible limit.

The objective is to model the flow and mixing behavior in a 2D mixing chamber of size `Lx = 3`, `Ly = 2` with three inlet flows and two rotating absorbing disks. The simulation solves:

- The non-dimensional continuity equation  
- The 2D incompressible Navier-Stokes equations  
- A convection-diffusion transport equation for a passive scalar representing species concentration (mass fraction `Y`)
