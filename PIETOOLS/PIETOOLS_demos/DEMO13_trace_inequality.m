% DEMO13_trace_inequality.m
%
% This document illustrates how a trace (boundary value) inequality can be
% verified, and its constant found, using PIETOOLS. The trace constant is
% the smallest value c such that
%   x(0)^2 + x(1)^2 ≤ c (||x_s||^2 + ||x||^2)    for all x in H1[0,1] = {x: x_s in L_2[0,1]},
% where there are no boundary conditions on x. We can pose this as an
% optimization program
%   min_{c}     c,
%      s.t.     x(0)^2 + x(1)^2 − c(<x_s,x_s> + <x,x>) <= 0
%               x \in H1.
% To solve, we declare the left boundary value as an input xa = x(0), so that
% the boundary value problem
%               x in W^1[0,1],     x(0) = xa,
% is solved by the free state xf = x_s in L2[0,1] and the free scalar w:
%               x(s) = xa + int_0^s x_f(theta) d theta.
% We find the operators such that
%               x(s) = H1op*[xa; xf] = xa + int_a^s x_s(theta) d theta ,
% and the operators that evaluate x at the boundary,
%               x(0) = E0op*[xa; xf] = xa,      
%               x(1) = E1op*[xa; xf] = xa + int_0^1 x_f(theta) d theta.
% Given these operators, we can pose the trace inequality as an LPI
%   min_{gam}   gam,
%      s.t.     gam*(H0op'*H0op + H1op'*H1op) - E0op'*E0op - E1op'*E1op >= 0,
% so that x(0)^2 + x(1)^2 <= gam(||x||^2 + ||x_s||^2) for all x in H1.
%
% Upper bound for c is max{4*lam, 4/lam} (lam = b - a), 
% trace theorem of Barraza Martinez et al., "On trace theorems and Poincare
% inequality for one-dimensional Sobolev spaces", arXiv:2108.13926 (2.5).
%
% The finite-element method gives bound 2.164 for [a, b] = [0, 1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - DEMO13
%
% Copyright (C)2026  PIETOOLS Team
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
% Initial coding, AT, 11/10/2026: modeled on DEMO3_poincare_inequality.

clc; clear; clear stateNameGenerator;
echo on

% =============================================
% === Declare the operators of interest

% % Declare the BVP
% %     x in W^1[a,b],   x(a) = w,   w in R free
a = 0;  b = 1;
pvar t s
x = pde_var(s,[a,b],1);     % x in W^1[a,b]
xa = pde_var('in',1);        % xa
BVP = subs(x,s,a)==xa;

% % Parameterize all x in H1:  x = Top*[w; xf],  xf = x_s
Top = solve_PIETOOLS_BVP(BVP,'all');
H0op = [Top.w{1}{1}, Top.x{1}{1}];  % (H0op*[w;xf]) = x
H1op = [Top.w{1}{2}, Top.x{1}{2}];  % (H1op*[w;xf]) = x_{s}
E0op = subs(H0op,s,a);              % (E0op*[w;xf]) = x(a)
E1op = subs(H0op,s,b);              % (E1op*[w;xf]) = x(b)


% =============================================
% === Declare the LPI

% % Initialize LPI program
prob = lpiprogram(H0op.vars,H0op.dom);

% % Declare decision variables:
% %   gam \in \R
[prob,gam] = lpidecvar(prob,'gam');     % scalar decision variable

% % Set inequality constraints:
% %   gam*(H0op'*H0op + H1op'*H1op) - E0op'*E0op - E1op'*E1op >= 0
prob = lpi_ineq_sop(prob,gam*(H0op'*H0op + H1op'*H1op) - E0op'*E0op - E1op'*E1op,struct('dD',2,'dw',2));

% % Set objective function:
% %   min gam
prob = lpisetobj(prob, gam);

% % Solve LPI and retrieve solution
prob = lpisolve(prob);
trace_constant = double(lpigetsol_sop(prob,gam));


echo off

fprintf(['\n If successful, ',num2str(trace_constant,4),' is an upper bound on c in x(0)^2+x(1)^2 <= c(||x||^2+||x_s||^2) for x in H1[0,1].\n'])
fprintf([' The optimal value is 2.164 (finite elements); the bound of Barraza Martinez et al. is 4.\n']);
