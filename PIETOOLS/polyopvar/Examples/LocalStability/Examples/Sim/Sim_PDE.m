
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Sim_PDE.m
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% CR, 10/07/2026: Initial coding.


%% Simulate PDE.

clear all; clear stateNameGenerator; close all; clc;

% Choose PDE and maximum simulation time.
PDE_name = "Heat"; % "Fisher", "Heat", "Burgers".
T = 2.0;

if PDE_name == "Fisher"
    % Define domain, parameters, and BCs of PDE being simulated.
    dom            = [0,1];
    L              = dom(2);
    alpha          = 5;
    beta           = -1;
    BC.left.type   = 'Dirichlet';
    BC.left.value  = 0;
    BC.right.type  = 'Dirichlet';
    BC.right.value = 0;

    % Define initial condition (which satisfies the BCs).
    %    u(x)         = q*sin(pi*x/L))
    %   ||u||_{Linf}  = abs(q)
    %   ||u||_{L2}    = abs(q) * sqrt(L/2)
    %   ||u_x||_{L2}  = abs(q) * ( pi/sqrt(2*L) )
    %   ||u||_{H1}    = abs(q) * sqrt( (L^2 + pi^2) / (2*L) )

    % q=5.7 is largest stable value for this initial condition.
    q  = 5.7;
    u0 = @(x) q * sin( (pi*x) / L );
    
    % Norm checks.
    Uinit_Linf = abs(q);
    Uinit_L2   = abs(q) * sqrt(L/2);
    Uxinit_L2  = abs(q) * ( pi/sqrt(2*L) );
    Uinit_H1   = abs(q) * sqrt( (L^2 + pi^2) / (2*L) );

elseif PDE_name == "Heat"
    % Define domain, parameters, and BCs of PDE being simulated.
    dom            = [0,1];
    L              = dom(2);
    a0             = 1;
    k              = 2*a0*pi^2 / L^2;
    BC.left.type   = 'Neumann';
    BC.left.value  = 0;
    BC.right.type  = 'Dirichlet';
    BC.right.value = 0;

    % Define initial condition (which satisfies the BCs).
    %    u(x)         = q * cos( (pi*x) / (2*L) )
    %   ||u||_{Linf}  = abs(q)
    %   ||u||_{L2}    = abs(q) * sqrt(L/2)
    %   ||u_x||_{L2}  = abs(q) * ( pi / (2*sqrt(2*L)) )
    %   ||u||_{H1}    = abs(q) * sqrt( (4*L^2 + pi^2) / (8*L) )
    
    % q=15.9 is largest stable value for this initial condition.
    % q=< a0*pi^2/L^2 can be verified by energy functional.
    q  = 15.9;
    u0 = @(x) q * cos( (pi*x) / (2*L) );
    
    % Norm checks.
    Uinit_Linf = abs(q);
    Uinit_L2   = abs(q) * sqrt(L/2);
    Uxinit_L2  = abs(q) * ( pi / (2*sqrt(2*L)) );
    Uinit_H1   = abs(q) * sqrt( (4*L^2 + pi^2) / (8*L) );

elseif PDE_name == "Burgers"
    % Define domain, parameters, and BCs of PDE being simulated.
    dom            = [0,1];
    L              = dom(2);
    v              = 1.0;
    r              = pi^2-0.1;
    BC.left.type   = 'Dirichlet';
    BC.left.value  = 0;
    BC.right.type  = 'Dirichlet';
    BC.right.value = 0;

    % Define initial condition (which satisfies the BCs).
    %   u(x)          = q*sin(pi*x/L))
    %   ||u||_{Linf}  = abs(q)
    %   ||u||_{L2}    = abs(q) * sqrt(L/2)
    %   ||u_x||_{L2}  = abs(q) * ( pi/sqrt(2*L) )
    %   ||u||_{H1}    = abs(q) * sqrt( (L^2 + pi^2) / (2*L) )
    
    % q=500 stable when r=pi^2-0.1.
    q  = 500;
    u0 = @(x) q * sin( (pi*x) / L );
    
    % Norm checks.
    Uinit_Linf = abs(q);
    Uinit_L2   = abs(q) * sqrt(L/2);
    Uxinit_L2  = abs(q) * ( pi/sqrt(2*L) );
    Uinit_H1   = abs(q) * sqrt( (L^2 + pi^2) / (2*L) );
else
    error("Sim_PDE: PDE_name set as unknown value.");
end

% Space and time grid.
dx = 0.005;
dt = 0.45*dx^2;

% Number of spatial grid points and time steps.
Nx = ceil((dom(2)-dom(1))/dx);
Nt = ceil(T/dt);

% Run the simulation.
if PDE_name == "Fisher"
    [x,t,U] = Fisher1D(alpha,beta,dom,T,Nx,Nt,u0);
elseif PDE_name == "Heat"
    [dx,x,t,U] = Heat1D(a0,k,dom,T,Nx,Nt,u0);
    BC.dx      = dx;
elseif PDE_name == "Burgers"
    [x,t,U] = Burgers1D(v,r,dom,T,Nx,Nt,u0);
end

x = x(:); % ensures x is column vector.

% Testing BCs hold.
[pass, err, tol] = BCs_test(U, BC);
if ~pass
    error('BCs_test:BoundaryConditionViolation', ...
        ['Boundary conditions are not satisfied. ' ...
         'Maximum left BC error = %.3e, ' ...
         'maximum right BC error = %.3e, ' ...
         'tolerance = %.3e.'], ...
         err.left, err.right, tol);
end

% Compute norms of U(t,x).
U_Linf  = zeros(1,Nt);
U_L2    = zeros(1,Nt); 
Ux_L2   = zeros(1,Nt);
U_H1    = zeros(1,Nt);
for i=1:Nt
    U_t       = U(:,i); % solution at that time.
    U_diff    = diff(U_t)./diff(x); % associated derivative.
    x_mid     = (x(1:end-1) + x(2:end))/2;

    U_Linf(i) = max(abs(U_t));
    U_L2(i)   = trapz(x, U_t.^2);
    Ux_L2(i)  = trapz(x_mid, U_diff.^2);
    U_H1(i)   = sqrt(U_L2(i) + Ux_L2(i));
    U_L2(i)   = sqrt(U_L2(i));
    Ux_L2(i)  = sqrt(Ux_L2(i));
end

%% Plots.

% Testing initial conditions match.
figure(1);
u_exact = u0(x);
plot(x, u_exact, 'k-', 'LineWidth', 2);
hold on;
plot(x, U(:,1), 'ro--');
xlabel('x');
legend('u_0(x)', 'U(:,1)', 'Interpreter', 'latex');

% 3D plot of solution.
figure(2);
[X,T] = meshgrid(x,t);
surf(T,X,U')
shading interp
xlabel('t');
ylabel('x');
zlabel('U(x,t)');
colormap(jet)

% Plot L_infty norm of U versus time.
figure(3);
plot(t,U_Linf(:),'r','LineWidth',2)
grid on
xlabel('t')
ylabel('$\|u\|_{L_{\infty}}$', 'Interpreter', 'latex');

% Plot L2 norm of U versus time.
figure(4);
plot(t,U_L2(:),'r','LineWidth',2)
grid on
xlabel('t')
ylabel('$\|u\|_{L_{2}}$', 'Interpreter', 'latex');

% Plot L2 norm of U_x versus time.
figure(5);
plot(t,Ux_L2(:),'r','LineWidth',2) 
grid on
xlabel('t')
ylabel('$\|u_{x}\|_{L_{2}}$', 'Interpreter', 'latex');

% Plot H1 norm of U versus time.
figure(6);
plot(t,U_H1(:),'r','LineWidth',2) 
grid on
xlabel('t')
ylabel('$\|u\|_{H_{1}}$', 'Interpreter', 'latex');