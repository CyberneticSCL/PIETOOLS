function [x,t,U] = Fisher1D(alpha,beta,dom,T,Nx,Nt,u0)
% [X,T,U] = FISHER1D(ALPHA,BETA,DOM,T,NX,NT,U0) simulates the Fisher
% equation
%   u_t = u_xx + alpha*u - beta*u^2,
% with u the population density, using an explicit Euler finite-difference
% scheme with homogeneous Dirichlet boundary conditions.
%
% INPUTS
% - alpha:  scalar, per-capita growth rate (alpha>0: the population
%           grows, alpha<0: the population decreases);
% - beta:   scalar, beta = alpha/K, with K the maximal sustainable
%           population density;
% - dom:    1 x 2 array [a,b] specifying the spatial interval;
% - T:      scalar, final time;
% - Nx:     integer, number of spatial grid points;
% - Nt:     integer, number of time steps;
% - u0:     function handle for the initial condition, u0(x), which is
%           assumed to satisfy the boundary conditions;
%
% OUTPUTS
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

% Time stepping (explicit Euler) with spatial loop defined to satisfy 
% Dirichlet BCs.
for n = 1:Nt-1
    for i = 2:Nx-1
        uxx      = (U(i+1,n) - 2*U(i,n) + U(i-1,n))/dx^2;
        reaction = alpha*U(i,n) - beta*U(i,n)^2;
        U(i,n+1) = U(i,n) + dt*(uxx + reaction);
    end
    
end
end