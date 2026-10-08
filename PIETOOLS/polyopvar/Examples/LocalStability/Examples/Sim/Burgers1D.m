function [x,t,U] = Burgers1D(v,r,dom,T,Nx,Nt,u0)
% Simulates:
% u_t = v*u_ss + r*u - u*u_s;
% u: velocity.
% v: viscosity (>=0).
% r: linear reaction (>0 --> growth; <0 --> decay).
%
% Inputs:
% v, r           - Burgers 1D parameters
% dom            - spatial interval
% T              - final time
% Nx             - number of spatial grid points
% Nt             - number of time steps
% u0             - function handle for initial condition, u0(x)
%
% Outputs:
% x  - spatial grid
% t  - time grid
% U  - solution matrix (Nx x Nt)

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
        ux  = (U(i+1,n) - U(i-1,n))/(2*dx);
        uxx = (U(i+1,n) - 2*U(i,n) + U(i-1,n))/dx^2;
        U(i,n+1) = U(i,n) + dt*(v*uxx + r*U(i,n) - U(i,n)*ux);
    end
    
end
end