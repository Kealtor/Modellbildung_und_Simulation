function MuS_FDM_2d_en
%
% Modeling and Simulation
% Exercise to chapter: Models with distributed parameters
%
% Balazs Pritz pritz@kit.edu
%
% 2D test case: lid-driven cavity flow
%
% Discretisation scheme: finite difference discretization methods
%                        central difference scheme
%
% Source: Griebel, Michael: Numerische Simulation in der Strömungsmechanik,
% Vieweg, 1995

clear all

% Time step
dt=0.1;

% Maximum number of time steps
% (Simulation time = timestep*maxnt)
% (It is possible to terminate the calculation through an other variable than maxnt,
% you can compute the change in the solutions of two subsequent time steps
% and use it as a criteria to terminate the calculation.)
maxnt=400;

% Velocity of the lid (wall at y=1)
u0=10.0;

% The convective terms are discretized with a combination of central differences
% and a modified upwind scheme.
% The blending parameter is gamma (gamma is here not the diffusion coefficient!).
% gamma=1 -> pure upwind diff. scheme
% gamma=0 -> pure central diff. scheme
gamma=0.8;


% Spatial resolution
% Domain of interest: x=[0,1], y=[0,1]
% Indicies: the walls are at i=2, i=nx+1 and at j=1, j=ny+1
% the first grid nodes in the flow are at i=3, i=nx, j=3, j=ny
% Dummy-nodes (1 and nx+2) willl be used to define boundary conditions for the velocity
%
% We use here an equidistant grid.
% If you want to use local refinement, you have to define vectors for x(i) and for y(i)
% and in the equations use x(i,j) und y(i,j) and replace dx and dy by
% x(i)-x(i-1), x(i+1)-x(i) or x(i+1)-x(i-1)
% (obviously dy by y()-y()).
%

% Number of grid nodes
nx=21;
ny=21;


nxm1=nx-1;
nym1=ny-1;
nxp1=nx+1;
nyp1=ny+1;
nxp2=nx+2;
nyp2=ny+2;


dx=1/(nxm1);
dy=1/(nym1);

% Material properties:
% Density
%  Air
rho=1.25;
%  Water
%rho=998;
% Dynamic viskosity
%  Air
mue=0.00001808;
%  Water
%mue=0.001895;

nue=mue/rho;

% Initialization
u=zeros(nxp2,nyp2);
v=zeros(nxp2,nyp2);
p=zeros(nxp2,nyp2);

% Coefficients for the derivatives at the wall
ew=ones(nxp1,1);
eo=ones(nxp1,1);
es=ones(nyp1,1);
en=ones(nyp1,1);
ew(2)=0;
eo(nxp1)=0;
es(2)=0;
en(nyp1)=0;

rhs=zeros(nxp2,nyp2);

% Time loop
for nt=1:maxnt
    
    % Boundary conditions for the velocity
    for i=2:nxp1
        u(i,1)=-u(i,2);
        u(i,nyp2)=2.*u0-u(i,nxp1);
        v(i,1)=0;
        v(i,nyp1)=0;
    end
    for j=2:nyp1
        u(1,j)=0;
        u(nxp1,j)=0;
        v(1,j)=-v(2,j);
        v(nxp2,j)=-v(nxp1,j);
    end

    % convective und viscous terms
    F=zeros(nxp1,nyp1);
    G=zeros(nxp1,nyp1);
    for j=2:nyp1
        for i=2:nx
            %
            duudx=0.25*( (u(i,j)+u(i+1,j))^2 - (u(i-1,j)+u(i,j))^2 + gamma*( abs(u(i,j)+u(i+1,j))*(u(i,j)-u(i+1,j)) - abs(u(i-1,j)+u(i,j))*(u(i-1,j)-u(i,j)) ) )/dx;
            duvdy=0.25*( (v(i,j)+v(i+1,j))*(u(i,j)+u(i,j+1)) - (v(i,j-1)+v(i+1,j-1))*(u(i,j-1)+u(i,j)) + gamma*( abs(v(i,j)+v(i+1,j))*(u(i,j)-u(i,j+1)) - abs(v(i,j-1)+v(i+1,j-1))*(u(i,j-1)-u(i,j)) ) )/dy;
            d2udx2=(u(i+1,j)-2.*u(i,j)+u(i-1,j))/dx^2;
            d2udy2=(u(i,j+1)-2.*u(i,j)+u(i,j-1))/dy^2;
            %
            % term for the time diskretization
            F(i,j)=dt*(nue*(d2udx2+d2udy2) - duudx - duvdy) + u(i,j);
        end
    end
    for j=2:ny
        for i=2:nxp1
            %
            duvdx=0.25*( (u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j)) - (u(i-1,j)+u(i-1,j+1))*(v(i-1,j)+v(i,j)) + gamma*( abs(u(i,j)+u(i,j+1))*(v(i,j)-v(i+1,j)) - abs(u(i-1,j)+u(i-1,j+1))*(v(i-1,j)-v(i,j)) ) )/dx;
            dvvdy=0.25*( (v(i,j)+v(i,j+1))^2 - (v(i,j-1)+v(i,j))^2 + gamma*( abs(v(i,j)+v(i,j+1))*(v(i,j)-v(i,j+1)) - abs(v(i,j-1)+v(i,j))*(v(i,j-1)-v(i,j)) ) )/dy;
            d2vdx2=(v(i+1,j)-2.*v(i,j)+v(i-1,j))/dx^2;
            d2vdy2=(v(i,j+1)-2.*v(i,j)+v(i,j-1))/dy^2;
            %
            % term for the time diskretization
            G(i,j)=dt*(nue*(d2vdx2+d2vdy2) - duvdx - dvvdy) + v(i,j);
        end
    end

    
    % Poisson equation for calculating the pressure field
    
    % The relaxation factor should be between 0 and 2
    % usually 1.7 is used, if it is 1, we get the Gauß-Seidel method
    omega=1.7;
    
    % maximum number of iterations solving the Poisson equation
    pitmax=100;
    
    for pit=1:pitmax
        p_old=p;

        %err=zeros(nxp1,nyp1);

        % Update boundary conditions
        for i=2:nxp1
            p(i,1)=p(i,2);
            p(i,nyp2)=p(i,nyp1);
            G(i,1)=v(i,1);
            G(i,nyp1)=v(i,nyp1);
        end
        for j=2:nyp1
            p(1,j)=p(2,j);
            p(nxp2,j)=p(nxp1,j);
            F(1,j)=u(1,j);
            F(nxp1,j)=u(nxp1,j);
        end

        for j=2:nyp1
            for i=2:nxp1
                rhs(i,j)=((F(i,j)-F(i-1,j))/dx + (G(i,j)-G(i,j-1))/dy)/dt;
                p(i,j)=(1.-omega)*p_old(i,j) + omega*((eo(i)*p_old(i+1,j)+ew(i)*p(i-1,j))/dx^2 + (en(j)*p_old(i,j+1)+es(j)*p(i,j-1))/dy^2 - rhs(i,j))/((eo(i)+ew(i))/dx^2 + (en(j)+es(j))/dy^2);
%                err(i,j)=(eo(i)*(p(i+1,j)-p(i,j)) - ew(i)*(p(i,j)-p(i-1,j)))/dx^2 + (en(j)*(p(i,j+1)-p(i,j)) - es(j)*(p(i,j)-p(i,j-1)))/dy^2 - rhs(i,j);
            end
        end
        % Convergence will be evaluated by calculating err or errit
%        errit=0;
%        for j=2:nyp1
%            for i=2:nxp1
%                errit=errit+err(i,j)^2;
%            end
%        end
%        errit=errit/(nxp1*nyp1);
%        % errit is the summed up error of the current iteration
%        % Here you could define a termination criterion
%        if errit < 
%        end
    end
    
    % New values for the velocity in the time will be calculated with an
    % explicte Euler scheme
    for j=2:nyp1
        for i=2:nx
            
            dpdx=(p(i+1,j)-p(i,j))/dx;
            
            % Term for the time discretization
            u(i,j)=F(i,j) - dt*dpdx/rho;
            
        end
    end
    for j=2:ny
        for i=2:nxp1
            
            dpdy=(p(i,j+1)-p(i,j))/dy;
            
            % Term for the time discretization
            v(i,j)=G(i,j) - dt*dpdy/rho;
            
        end
    end
    
    % If CFL becomes big, the time derivative could be unstable.
    % In such a situation you should reduce the time step size.
    dcfl=0;
    for j=2:ny
        for i=2:nxp1
            
            dcfl=max(2*nue*dt/(dx*dx) + max(u(i,j)/dx,v(i,j)/dy)*dt,dcfl);
            
        end
    end
     
    % Simulation time
    zeit=nt*dt;

    % Output
    % The solution matrices have to be transformed
    % Graphical output will be generated after each n-th iteration
    % (the generation of a plot takes computation time,
    % if you set n higher, your simulation runs faster).
    % Pressure will be plotted as contour, the velocities as vectors.
    n=10;
    if rem(nt,n)==0
        % Update the boundary conditions
        for i=2:nxp1
            u(i,1)=-u(i,2);
            u(i,nyp2)=2.*u0-u(i,nxp1);
            v(i,1)=0;
            v(i,nyp1)=0;
        end
        for j=2:nyp1
            u(1,j)=0;
            u(nxp1,j)=0;
            v(1,j)=-v(2,j);
            v(nxp2,j)=-v(nxp1,j);
        end
    
        pp=zeros(ny,nx);
        up=zeros(ny,nx);
        vp=zeros(ny,nx);
        for j=1:ny
            for i=1:nx
                pp(j,i)=p(i+1,j+1);
                up(j,i)=u(i+1,j+1);
                vp(j,i)=v(i+1,j+1);
            end
        end
        figure(1)
        [X,Y] = meshgrid(0:dx:1,0:dy:1);
        contourf(X,Y,pp)
        colorbar
        hold on
        quiver(X,Y,up,vp,'w')
        text(0,1.05,['DCFL= ',num2str(dcfl)])
        text(0.3,1.05,['Zeit= ',num2str(zeit)])
        hold off
        
    end
end
%End of time loop

end
% End of function

