% PIESOS_DEMO2_local_stability.m
%
% This document illustrates how PIESOS can be used to verify local
% stability of a nonlinear PDE on a ball in the L2 norm. Specifically, we
% consider the Fisher equation:
%   u_{t}(t,s) = u_{ss}(t,s) + r*u(t,s) + u(t,s)^2,     s in [0,1];
%       u(t,0) = u(t,1) = 0,
% which is the Fisher-KPP equation v_{t} = v_{ss} + r*v*(1-v) for v=-u/r.
% The zero solution is stable if and only if r <= pi^2 = 9.8696, but not
% globally: for r < pi^2 the PDE has an unstable steady state u*, with
% u*_{ss} + r*u* + u*^2 = 0, whose L2 norm ||u*|| shrinks to zero as r goes
% to pi^2. For r = 9, we have ||u*||=0.7238, so stability can be verified
% on a ball of at most R<0.7238
%
% We test stability on the ball ||u|| <= R, by constructing a Lyapunov
% functional V(u) = <u,P*u> + eppos*||u||^2, with P(s) >= 0 a polynomial
% multiplier, and enforcing the distributed sum of squares (SOS) conditions
%   V <= gam*||u||^2,                       (gam is minimized)
%   dV + lam*(R^2-||u||^2) <= 0,            (lam >= 0 is an SOS multiplier)
% with dV the derivative of V along solutions to the PDE. By the second
% condition, dV <= 0 whenever ||u|| <= R, so that solutions starting in the
% ball of radius rho = R*sqrt(eppos/gam) remain in the ball of radius R.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESOS_DEMO2
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

% Declare independent variables (time and space), the reaction rate r, and
% the radius R of the ball on which to test stability
pvar t s
dom = [0,1];
r = 9;
R = 0.72;       % for r=9, the largest stable ball has radius at most ||u*||=0.7238

% Declare the PDE state, and the PDE with its boundary conditions
x = pde_var(s,dom);
PDE = [diff(x,t)==diff(x,s,2)+r*x+x^2;
       subs(x,s,dom(1))==0;
       subs(x,s,dom(2))==0];

% Convert to a PIE: T*x_{t} = f(x), with u = T*x the PDE state
PIE = convert(PDE,'silent');
T = PIE.T;      f = PIE.f;      x = PIE.x;
u = T*x;

% Declare the ball ||u||<=R as the set where g(u) = R^2-<u,u> >= 0
g = R^2 - u'*u;


% =============================================
% === Declare the PIESOS program

% % Initialize the program, with the PIE state x and the upper bound gam on
% % the Lyapunov functional as variables, and minimize gam
dpvar gam
prog = piesos_program(x,gam);
prog = piesos_setobj(prog,gam);

% % Declare the Lyapunov functional V = <u,P*u> + eppos*<u,u>, with
% % P(s)=Z(s)'*Q*Z(s)>=0 a multiplier of degree pdeg, and compute the
% % derivative dV of V along solutions to the PIE
pdeg = 3;               % degree of the monomials in s used to parameterize P
eppos = 1;              % V >= eppos*<u,u>
Popts.degs = pdeg;
Popts.exclude = [0;1;1];    % multiplier operator only
Popts.pos = true;           % V >= 0
[prog,V,dV] = piesos_LF(prog,PIE,'pde',1,Popts);
V  = V + eppos*(u'*u);
dV = dV + 2*eppos*(u'*f);

% % Enforce the upper bound V <= gam*<u,u>
ineq_opts.degs = pdeg;      % degree of the monomials in s of the SOS terms
ineq_opts.psatz = 1;        % allow multipliers (s-a)*(b-s)
prog = piesos_ineq(prog,gam*(u'*u)-V,ineq_opts);

% % Declare an SOS multiplier lam of degree 2 in u, and enforce
% %   -dV - lam*g >= 0
% % so that dV <= 0 whenever g >= 0
lam_opts.degs = pdeg;
lam_opts.exclude = [1,0,0]';    % integral operators only
[prog,lam] = piesos_sosvar(prog,dmonomials(x,1),lam_opts);
ineq_opts.degs = pdeg+1;
prog = piesos_ineq(prog,-(dV+lam*g),ineq_opts);

% % Solve and retrieve the solution
% sol_opts.solver = 'sedumi';   % uncomment to declare SDP solver
sol_opts.simplify = true;       % simplify SDP before solving
prog_sol = piesos_solve(prog,sol_opts);

echo off

% Check if the PDE is stable on the ball ||u||<=R
is_pinf = prog_sol.solinfo.info.pinf;       % is primal feasible?
is_dinf = prog_sol.solinfo.info.dinf;       % is dual feasible?
is_nerr = prog_sol.solinfo.info.numerr;     % numerical errors?
feasrat = prog_sol.solinfo.info.feasratio;  % ratio should be close to 1
if is_pinf || is_dinf || is_nerr || abs(feasrat-1)>0.1
    fprintf('\n Stability on the ball of radius R=%g could not be verified.\n',R)
else
    gam_sol = double(piesos_getsol(prog_sol,gam));
    fprintf('\n Local stability of the PDE was verified for r=%g and R=%g.\n',r,R)
    fprintf(' Solutions starting in the ball of radius %g remain in the ball of radius %g.\n', ...
        R*sqrt(eppos/gam_sol),R)
end
fprintf('\n For r=9, a non-trivial equilibrium exists with ||u||_{L2}=0.7238.\n')