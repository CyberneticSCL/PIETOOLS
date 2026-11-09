% PIESOS_DEMO1_global_stability.m
%
% This document illustrates how PIESOS can be used to verify global
% stability of a nonlinear PDE, using the energy V(u) = <u,u> = ||u||^2 as
% Lyapunov functional. Specifically, we consider Burgers' equation with an
% added reaction term:
%   u_{t}(t,s) = u_{ss}(t,s) + r*u(t,s) - u(t,s)*u_{s}(t,s),    s in [0,1];
%       u(t,0) = u(t,1) = 0.
% Along solutions to this PDE, the energy satisfies
%   dV = 2*<u,u_{ss}> + 2*r*||u||^2 - 2*int_0^1 u(s)^2*u_{s}(s) ds
%      = -2*(||u_{s}||^2 - r*||u||^2),
% since the nonlinear term integrates to zero: 
%   int_{0}^{1} u^2*u_{s} ds = [u^3/3]_0^1.
% The PDE is therefore globally stable if 
%   ||u_{s}||^2 - r*||u||^2 >= 0 
% for all u, which holds if and only if r <= pi^2 = 9.8696 (by Poincare).
% For r > pi^2, the PDE is unstable, as its linearization is unstable.
%  
% The test below enforces this positivity condition as a distributed sum
% of squares (SOS) constraint, -dV = W >= 0, using PIESOS. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESOS_DEMO1
%
% Copyright (C) 2026  PIETOOLS Team
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
% DJ, 11/07/2026: Initial coding;

clc; clear; clear stateNameGenerator
echo on

% =============================================
% === Declare the nonlinear PDE

% Declare independent variables (time and space), and the reaction rate
pvar t s
dom = [0,1];
r = pi^2 - 0.001;       % the PDE is stable if and only if r <= pi^2

% Declare the PDE state, and the PDE with its boundary conditions
x = pde_var(s,dom);
PDE = [diff(x,t)==diff(x,s,2)+r*x-x*diff(x,s);
       subs(x,s,dom(1))==0;
       subs(x,s,dom(2))==0];

% Convert to a PIE: T*x_{t} = f(x), with u = T*x the PDE state
PIE = convert(PDE,'silent');
T = PIE.T;      f = PIE.f;      x = PIE.x;


% =============================================
% === Declare the Lyapunov functional and its derivative

% Declare the energy V = <u,u>, with u = T*x the PDE state, and compute its
% derivative along solutions to the PIE: dV = 2*<u,f(x)>
u = T*x;
V  = u'*u;
dV = 2*(u'*f);


% =============================================
% === Declare the PIESOS program

% % Initialize the program, with the PIE state x as variable
prog = piesos_program(x);

% % Enforce -dV >= 0 on all u, as a distributed SOS constraint -dV = W:
% %   W = <Zop*x, Q*Zop*x>,   Q >= 0
% %   degs:   degree of the monomials in s used to parameterize Q
% %   psatz:  allow multipliers (s-a)*(b-s) in the parameterization of Q
ineq_opts.degs = 5;
ineq_opts.psatz = 1;
[prog,W] = piesos_ineq(prog,-dV,ineq_opts);

% % Solve and retrieve the solution
% sol_opts.solver = 'sedumi';   % uncomment to declare SDP solver
sol_opts.simplify = true;       % simplify SDP before solving
prog_sol = piesos_solve(prog,sol_opts);

echo off

% Check if the PDE is globally stable: dV = -W <= 0 for all u
is_pinf = prog_sol.solinfo.info.pinf;       % is primal feasible?
is_dinf = prog_sol.solinfo.info.dinf;       % is dual feasible?
is_nerr = prog_sol.solinfo.info.numerr;     % numerical errors?
feasrat = prog_sol.solinfo.info.feasratio;  % ratio should be close to 1
if is_pinf || is_dinf || is_nerr || abs(feasrat-1)>0.1
    fprintf('\n Global stability could not be verified for r=%g.\n',r)
else
    fprintf('\n Global stability of the PDE was verified for r=%g.\n',r)
end
fprintf(' The PDE is stable if and only if r<=pi^2=%g.\n',pi^2);
