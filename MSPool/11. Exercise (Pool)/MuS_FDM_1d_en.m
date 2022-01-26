function MuS_FDM_1d_en
%
% Modeling and Simulation
% Exercise to chapter: Models with distributed parameters
%
% Dr.-Ing. Balazs Pritz, pritz@kit.edu
%
% 1D test case: convection/diffusion equation with Dirichlet boundary conditions
%
% Discretisation scheme: finite difference discretization methods
% 
% Source: Joel H. Ferziger, Milovan Peric:
% Computational Methods for Fluid Dynamics, Springer, 2002

clear all

% time step
dt=0.0000193;

% Maximum number of time steps
% (Simulation time = timestep*maxnt)
% (It is possible to terminate the calculation through an other variable than maxnt,
% you can compute the change in the solutions of two subsequent time steps
% and use it as a criteria to terminate the calculation.)
maxnt=300;

% Alternativly you can define the simulation time and calculate maxnt based on the time step.

% Convection velocity
u0=4;

% Spatial resolution
% Domain of interest: x=[0,1]
%
% We use here an equidistant grid.
% If you want to use local refinement, you have to define a vector for x(i)
% and in the equations replace dx by
% x(i)-x(i-1), x(i+1)-x(i) or x(i+1)-x(i-1).
%
% Number of grid nodes
nx=500;

nxm1=nx-1;

% Cell size
dx=1/(nxm1);

% Material properties:
% Density
%  Air
rho=1;
%  Water
%rho=998;

% Diffusion coefficient for phi
gamma=0.1;


% Initialization
phi_uds=zeros(nx,1);
phi_cds=zeros(nx,1);

% Boundary conditions for phi
phi_uds(1)=0;
phi_uds(nx)=1;

phi_cds(1)=0;
phi_cds(nx)=1;


% Peclet-number
% If you want to set a particular Pe, you can define freely rho, u or gamma
peclet=rho*u0/gamma;


% The exact solution (for this we use a better spatial resolution)
xel=zeros(100,1);
el=zeros(100,1);
for i=1:100
    xel(i)=0.01*i-0.01;
    el(i)=(exp(xel(i)*peclet) - 1.)/(exp(peclet) - 1.);
end

% To evaluate the stability DCFL will be calculated.
dcfl=0;
for i=2:nx
    %
    dcfl=max(2*gamma*dt/(dx*dx*rho) + u0*dt/dx,dcfl);
    %
end


% Spatial position of the grid nodes
x=zeros(nx,1);
for i=1:nx
    x(i)=dx*i-dx;
end


% Initialization of vectores
rhs_uds=zeros(nx,1);
rhs_cds=zeros(nx,1);


% Time loop
for nt=1:maxnt
    
    for i=2:nxm1
        % Convective term with upwind difference scheme
        konv_uds=rho*u0*((phi_uds(i)-phi_uds(i-1))/dx);  % >>> uds dphi/dx comes here for phi_uds <<<;
        
        % convective term with central difference scheme
        konv_cds=rho*u0*((phi_cds(i+1)-phi_cds(i-1))/(2*dx)); % >>> cds dphi/dx comes here for phi_cds<<<;
        
        % Diffusive term with CDS
        diff_uds=gamma*(phi_uds(i+1)-2*phi_uds(i)+phi_uds(i-1))/dx^2 ; % >>> cds d2phi/dx2 comes here for phi_uds <<<;
        diff_cds=gamma*(phi_uds(i+1)-2*phi_uds(i)+phi_uds(i-1))/dx^2 ; % >>> cds d2phi/dx2 comes here for phi_cds <<<;

        rhs_uds(i)=diff_uds-konv_uds;
        rhs_cds(i)=diff_cds-konv_cds;

    end
    
    % Calculation of new values of phi in the time with an
    % explicite Euler scheme
    for i=2:nxm1
        phi_uds(i)= phi_uds(i)+rhs_uds(i)*dt; % >>> expicit Euler comes here for phi_uds <<<;
        phi_cds(i)= phi_cds(i)+rhs_cds(i)*dt; % >>> expicit Euler comes here for phi_cds <<<;
    end

    % Simulation time
    zeit=nt*dt;

    % Graphical output of the solution
    % after each n-th iteration
    % (the generation of a plot takes computation time,
    % if you set n higher, your simulation runs faster)
    n=20;
    if rem(nt,n)==0
        figure(1)
        plot(xel,el,'-b',x,phi_uds,'-ro',x,phi_cds,'-kx');
        legend('exact','uds','cds','Location','Eastoutside');
        axis([0 1 -0.2 1])
        text(0.1,0.9,['Pe: ',num2str(peclet)])
        text(0.1,0.85,['DCFL: ',num2str(dcfl)])
        text(0.1,0.8,['Time= ',num2str(zeit)])
        text(0.1,0.75,['Timestep= ',num2str(dt)])
        text(0.1,0.7,['Number of Grids = ',num2str(nx)])
       % text(0.1,0.65,['Convection Velocity = ',num2str(u0)])
        %drawnow % in new version of matlab this command must be used
    end
   
end
%End of time loop

end
% End of function

