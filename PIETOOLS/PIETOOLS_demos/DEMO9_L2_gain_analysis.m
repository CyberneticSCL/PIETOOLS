% DEMO9_L2_gain_analysis.m
% See Chapter 11.9 of the manual for a description.
%
% This document illustrates how an upper bound on the L2-gain of a 2D PDE 
% with inputs and outputs can be computed using PIETOOLS
% Specifically, we consider the following system:
% PDE:      d/dt x(t,s1,s2) = (d^2/ds1^2) x(t,s) +(d^2/ds2^2) x(t,s) + lam*x(t,s1,s2) + w(t),   s1,s2 in [0,1];
% Outputs:             z(t) = int_{0}^{1} int_{0}^{1} x(t,s1,s2) ds1 ds2;
% BCs:                    0 = x(t,0,s2) = x(t,1,s2);
%                         0 = x(t,s1,0) = d/ds2 x(t,s1,1);
% (unstable for lam > 5/4 pi^2)
% First we convert the above PDE to a PIE of the form:
%   [T \dot{v}](t,s) = [A v](t,s) + [B w](t,s);
%               z(t) = [C v](t)  + [D w](t);
% Then, we find gam such that any solution (w,z) satisfies
%    ||z||_{L2[0,infty)}/||w||_{L2[0,infty)} <= gam
% by solving the LPI
%   min_{gam,P} gam
%   s.t.          P>=0,
%       [-gam*I,   -D',     -P*B'*T          ] <=0
%       [-D,       -gam*I,  C                ]
%       [-T'*P*B,  C',      (P*A)'*T+T'*(P*A)]
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - DEMO8
%
% Copyright (C)2024  PIETOOLS Team
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
% DJ, 11/22/2024: Initial coding;
% DJ, 12/16/2024: Simplify demo;
% DJ, 12/26/2024: Add simulation;
% DJ, 12/28/2024: Change simulation conditions;

clc; clear; close all; clear stateNameGenerator
echo on

% =============================================
% === Declare the system of interest

% % Declare system as PDE
% Declare independent variables (time and space)
pvar s1 s2 t
% Declare state, input, and output variables
a = 0;      b = 1;
c = 0;      d = 1;
x = pde_var('state',1,[s1;s2],[a,b;c,d]);
w = pde_var('in',1);
z = pde_var('out',1);
% Declare the sytem equations
lam = 5;
PDE = [diff(x,t) == diff(x,s1,2) +diff(x,s2,2) + lam*x + w;
       z == int(x,[s1;s2],[a,b;c,d]);
       subs(x,s1,a)==0;
       subs(x,s1,b)==0;
       subs(x,s2,c)==0;
       subs(diff(x,s2),s2,d)==0];
PDE = initialize(PDE);
display_PDE(PDE);

% % Convert PDE to PIE
PIE = convert(PDE);
T = PIE.T;
A = PIE.A;     C = PIE.C1;
B = PIE.B1;    D = PIE.D11;


% =============================================
% === Simulate the system

% % Declare initial values and disturbance
syms st sx sy real
uinput.ic.PDE = 0;
uinput.w = 20*sin(pi*st)*exp(-st/2);

% % Set options for discretization and simulation
opts.plot = 'yes';  % plot the final solution
opts.N = 8;         % Expand using 8x8 Chebyshev polynomials
opts.tf = 10;       % Simulate up to t = 10
opts.dt = 1e-2;     % Use time step of 10^-2
ndiff = [0,0,1];    % The state involves 1 second order differentiable state variable wrt sx and sy

% % Perform the actual simulation
[solution,grid] = PIESIM(PDE,opts,uinput,ndiff);
tval = solution.timedep.dtime;
x = reshape(solution.timedep.pde{2}(:,:,1,:),opts.N+1,opts.N+1,[]);
z = solution.timedep.regulated(1,:);
w = double(subs(uinput.w,st,tval));


% =============================================
% === Declare the LPI

% % Initialize LPI program
prog = lpiprogram([s1;s2],[a,b;c,d]);

% % Declare decision variables:
% %   gam \in \R,     P:L2-->L2,    Z:\R-->L2
% Scalar decision variable
[prog,gam] = lpidecvar(prog,'gam');
% Positive operator variable P>=0
[prog,P] = poslpivar(prog,T.dim);

% % Set inequality constraints:
% %   Q <= 0
Q = [-gam,        D',      (P*B)'*T;
     D,           -gam,    C;
     T'*(P*B),    C',     (P*A)'*T+T'*(P*A)];
opts_Q.psatz = 2;
prog = lpi_ineq(prog,-Q,opts_Q);

% % Set objective function:
% %   min gam
prog = lpisetobj(prog, gam);

% % Solve and retrieve the solution
prog_sol = lpisolve(prog,opts);
% Extract solved value of decision variable
gam_val = double(lpigetsol(prog_sol,gam));

% % Alternatively, uncomment to use pre-defined executive
% [prog, ~, gam_val] = lpiscript(PIE,'l2gain','stripped');   


echo off

fprintf(['\n If successful, ',num2str(gam_val,7),' is an upper bound on the L2-gain from disturbance to output for the 2D reaction-diffusion equation.\n'])
fprintf([' The true value of the L2-gain for this system with lam=5 and on [0,1]^2 is known to be roughly 0.09397.\n']);