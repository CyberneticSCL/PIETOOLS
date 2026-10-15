function [dx,x,t,U] = Heat1D(a0,k,dom,T,Nx,Nt,u0)
% [DX,X,T,U] = HEAT1D(A0,K,DOM,T,NX,NT,U0) simulates the semilinear heat
% equation
%   u_t = a0*u_xx + u^2 - (k/b)*\int_0^b u dx,
% with u the temperature, using an explicit Euler finite-difference scheme
% with a Neumann boundary condition u_x(a)=0 on the left and a Dirichlet
% boundary condition u(b)=0 on the right.
%
% INPUTS
% - a0:     scalar, diffusion coefficient (a0>0);
% - k:      scalar, integral feedback gain;
% - dom:    1 x 2 array [a,b] specifying the spatial interval;
% - T:      scalar, final time;
% - Nx:     integer, number of spatial grid points;
% - Nt:     integer, number of time steps;
% - u0:     function handle for the initial condition, u0(x). The
%           boundary conditions are enforced on the initial condition;
%
% OUTPUTS
% - dx:     scalar, spatial grid spacing;
% - x:      1 x Nx array, spatial grid;
% - t:      1 x Nt array, time grid;
% - U:      Nx x Nt array, solution matrix, with U(i,n) the solution at
%           x(i) and t(n);

% Grids
a  = dom(1);
b  = dom(2);
dx = (b-a)/(Nx-1);
dt = T/Nt;

x = linspace(a,b,Nx);
t = linspace(0,T,Nt);

% Stability condition for explicit diffusion.
if dt/dx^2 > 0.5
    warning('Scheme may be unstable: dt/dx^2 > 0.5')
end

% Allocate solution matrix.
U = zeros(Nx,Nt);

% Initial condition - assumes BCs are satisfied.
U(:,1)   = u0(x)';

% Enforce BCs on initial condition.
U(1,1)  = U(2,1);   % Neumann: u_x(0) = 0
U(Nx,1) = 0;        % Dirichlet: u(b) = 0

% Time stepping (explicit Euler).
for n = 1:Nt-1
    feedback = (k/b)*trapz(dx,U(:,n));
    uxx = (U(3,n) - 2*U(2,n) + U(1,n))/dx^2;
    U(1,n+1) = U(1,n) + dt*(a0*uxx + U(1,n)^2 - feedback); % Left Neumann BC.

    for i = 2:Nx-1
        uxx      = (U(i+1,n) - 2*U(i,n) + U(i-1,n))/dx^2;
        U(i,n+1) = U(i,n) + dt*(a0*uxx + U(i,n)^2 - feedback);
    end  
end


end