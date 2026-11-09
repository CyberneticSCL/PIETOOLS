% PIESOS_DEMO3_energy_bound.m
%
% This document illustrates how PIESOS can be used to bound a nonconvex
% functional from below. Specifically, we bound the Ginzburg-Landau energy
%   E(u) = int_0^1 u_{s}^2 - r*u^2 + u^4/2 ds,       u(0) = u(1) = 0,
% from below, over all functions u. This is the energy of the Allen-Cahn
% equation u_{t} = u_{ss} + r*u - u^3, for which dE/dt = -2*||u_{t}||^2.
% No PDE is solved, and no Lyapunov functional is sought: we only enforce
%   E(u) - gam >= 0   for all u,
% as a distributed sum of squares (SOS) constraint, for the largest gam.
% For r <= pi^2, E >= 0 by the Poincare inequality, and E(0) = 0, so that
% the minimum is 0. For r > pi^2, the minimum is attained by the positive
% solution u* of u_{ss} + r*u - u^3 = 0, with minimum 
%   E(u*) = -||u*||_4^4/2.
% We compute this exact minimum from the boundary value problem, and 
% compare it to the bound from the SOS constraint, for which gam is found
% by bisection.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESOS_DEMO3
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
% DJ, 11/08/2026: Initial coding;

clc; clear; clear stateNameGenerator
echo on

% =============================================
% === Declare parameters

% Declare independent variables (time and space), and the reaction rate
pvar t s
dom = [0,1];
r = 12;

% Set bisection limits for the lower bound gam on E. Since u_{s}^2 >=
% pi^2*u^2 (Poincare) and u^4/2 - c*u^2 >= -c^2/2, we have E >= -(r-pi^2)^2/2,
% and E(0) = 0, so that gam <= 0
gam_min = -(r-pi^2)^2/2;        gam_max = 0;
gam = 0.5*(gam_min + gam_max);
n_iters = 10;


% =============================================
% === Declare the energy

% Declare the Allen-Cahn PDE, only to obtain the PIE state x = u_{ss} and
% the map u = T*x from x to the PDE state
u = pde_var(s,dom);
BVP = [diff(u,t)==diff(u,s,2)+r*u-u^3;
       subs(u,s,dom(1))==0;
       subs(u,s,dom(2))==0];
Top = solve_PIETOOLS_BVP(BVP);
x = polyopvar('x',s,dom);
u = Top.x*x;

% Declare the energy, using int u_{s}^2 ds = -<u,u_{ss}> = -<u,x>
u2 = u.^2;
E = -(u'*x) - r*(u'*u) + 0.5*(u2'*u2);


%%% Perform bisection on the value of gam
for iter = 1:n_iters
    fprintf(['\n',' --- Testing E >= ',num2str(gam),' ---\n'])
    % =============================================
    % === Declare the PIESOS program

    % % Initialize the program, with the PIE state x as variable
    prog = piesos_program(x);

    % % Enforce E - gam >= 0 as a distributed SOS constraint:
    % %   E - gam = <Zop*Z(x), Q*Zop*Z(x)>,    Q >= 0,   Z(x) = [1; x; x^2]
    % %   degs:   degree of the monomials in s used to parameterize Zop
    % %   psatz:  allow multipliers (s-a)*(b-s) in the parameterization of Q
    ineq_opts.degs = 4;
    ineq_opts.psatz = 1;
    [prog,W] = piesos_ineq(prog,E-gam,ineq_opts);

    % % Solve and retrieve the solution
    sol_opts.simplify = true;       % simplify SDP before solving
    prog_sol = piesos_solve(prog,sol_opts);

    % Check if the bound is verified. The ratio feasrat is not close to 1
    % for feasible programs here, but is far from 1 for infeasible ones
    is_pinf = prog_sol.solinfo.info.pinf;       % is primal feasible?
    is_dinf = prog_sol.solinfo.info.dinf;       % is dual feasible?
    is_nerr = prog_sol.solinfo.info.numerr;     % numerical errors?
    feasrat = prog_sol.solinfo.info.feasratio;
    if is_pinf || is_dinf || is_nerr || abs(feasrat-1)>0.5
        % The bound cannot be verified --> decrease gam
        gam_max = gam;
        gam = 0.5*(gam_min + gam_max);
    else
        % The bound is verified --> try larger gam
        gam_min = gam;
        gam = 0.5*(gam_min + gam_max);
    end
end
echo off

fprintf('\n E(u) >= %g was verified for all u.\n',gam_min)
% For r>0, the exact minimum E(u*) = -||u*||_4^4/2 can be computed from the
% boundary value problem u_{ss} + r*u - u^3 = 0, using a guess 
% u = a*sin(pi*s)
fprintf(' For r=12, the exact minimum of E is -1.5218.\n')