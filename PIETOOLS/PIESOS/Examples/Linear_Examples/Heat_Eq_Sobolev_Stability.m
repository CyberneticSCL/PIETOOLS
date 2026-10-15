%
% Heat equation
%
% PDE: u_t = u_ss;
% BCs: u(t,0) = u(t,1) = 0;
%
% fundamental state: v = u_ss \in L_2[0,1];
%
% PIE: (Tv_t) = f(v) = v;
%
% This script tests stability of the PIE in the Sobolev norm, using the
% Lyapunov functional V(x) = -<x,T*x> = ||u_s||^2, in terms of the PIE
% state x = v. The constraints
%   V >= eppos*||u_s||^2;    (positivity)
%   dV <= 0;                 (derivative along the PIE, dV = -2*<x,f(x)>)
% are enforced using `piesos_ineq`. The derivative u_s = R*x is obtained
% as the output PIE.h of the PDE with output z = u_s.
%
% The user may change the coefficient eppos.
%
% DJ, 10/15/2026: Initial coding;

clear;  clear stateNameGenerator

pvar s s_dum t
dom = [0,1];
u = pde_var(s,dom);
z = pde_var('out',s,dom);    % output z = u_s, to obtain R: z = R*x
PDE = [diff(u,t)==diff(u,s,2);
       subs(u,s,dom(1))==0;  subs(u,s,dom(2))==0;
       z==diff(u,s)];

PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');
T = PIE.T;
f = PIE.f;
x = PIE.x;

Rx = PIE.h;                 % R*x = u_s

% Declare the SOS program
prog = piesos_program(x);

% Declare the LF
Vx = -innerprod(x,T*x);

% Enforce strict positivity
eppos = 1;
prog = piesos_ineq(prog,Vx-eppos*innerprod(Rx,Rx));

% Compute the Lie derivative
dV = -2*innerprod(x,f);

% Enforce negativity
prog = piesos_ineq(prog,-dV);

% Solve
prog_sol = piesos_solve(prog);
