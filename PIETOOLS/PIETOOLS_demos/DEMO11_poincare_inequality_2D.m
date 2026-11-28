% DEMO11_poincare_inequality_2D.m
%
% This document illustrates how the Poincare constant of a 2-D domain can
% be found using PIETOOLS. The Poincare constant is the smallest c such that
%   ||x|| <= c ( ||x_{s1}||^2 + ||x_{s2}||^2 )^{1/2}
% for all x(s1,s2) on [a,b]x[c,d] that are periodic in s2 and vanish at
% s1 = a,b. We pose this as an optimization program
%   min_{gam}   gam,
%      s.t.     <x,x> - gam (<x_{s1},x_{s1}> + <x_{s2},x_{s2}>) <= 0,
% where now c = sqrt(gam). We suppose x is twice differentiable
% in each variable, declare the boundary conditions as a BVP, and use
% CONVERT to find the operator H0op such that x(t) = H0op*xf(t), for a fundamental
% state xf; then H1op = d/ds1 H0op, H2op = d/ds2 H0op give x_{s1}, x_{s2}.
% The Dirichlet conditions at s1 = a,b can be solved out, but the periodic
% conditions in s2 cannot (the BVP is not uniquely solvable). They remain a
% constraint K*xf(t) = 0 on the fundamental state, which we impose using
% Finsler's lemma with a decision-variable slack X:
%   min gam  s.t.  -( H0'H0 - gam (H1'H1 + H2'H2) + X K + K' X' ) >= 0.
% The exact value is (b-a)/pi.
%
% The script solves the LPI on [-1,1]x[0,pi] and prints the bound
% sqrt(gam) next to (b-a)/pi. Workspace: H0op, H1op, H2op, Kop, the solved
% program prob, poincare_constant. ECHO ON prints each line as it runs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - DEMO11
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding DB - 11/18/2026

% Clear the workspace; print each line as it runs
clc; clear; clear stateNameGenerator;
echo on

% =============================================
% === Declare the operators of interest

% % Declare the BVP
% %   x(a,s2) = x(b,s2) = 0,   x(s1,c) = x(s1,d),   x_{s2}(s1,c) = x_{s2}(s1,d)
a = -1;  b = 1;  c = 0;  d = pi;
dom = [a b; c d];
pvar s1 s2
x = pde_var([s1;s2],dom,[2;2]);     % x in W^2 in both variables
BVP = [subs(x,s1,a)==0;
       subs(x,s1,b)==0;
       subs(x,s2,c)==subs(x,s2,d);
       subs(diff(x,s2),s2,c)==subs(diff(x,s2),s2,d)];

% % Convert the BVP:  x(t) = H0op*xf(t),  K*xf(t) = 0
BVP_sol = convert(BVP);
H0op = BVP_sol.T;                   % (H0op*xf(t)) = x(t)
H1op = diff(H0op,s1);               % (H1op*xf(t)) = x_{s1}(t)
H2op = diff(H0op,s2);               % (H2op*xf(t)) = x_{s2}(t)
Kop = BVP_sol.K;                    % constraint K*xf(t) = 0


% =============================================
% === Declare the LPI

% % Initialize LPI program
prob = lpiprogram(H0op.vars,H0op.dom);

% % Declare decision variables:
% %   gam \in \R,   X: constraint space --> space of xf
[prob,gam] = lpidecvar(prob,'gam');     % scalar decision variable
% X maps the constraint space to the space of xf, with the output and input
% spaces of Kop'; monomial degree 1
[prob,X] = lpivar(prob,Kop',1);

% % Set inequality constraint:
% %   gam*(H1'H1 + H2'H2) - H0'H0 - X*K - K'*X' >= 0
prob = lpi_ineq(prob,gam*(H1op'*H1op+H2op'*H2op)-H0op'*H0op-X*Kop-Kop'*X');

% % Set objective function:
% %   min gam
prob = lpisetobj(prob, gam);

% % Solve LPI and retrieve solution
prob = lpisolve(prob);
poincare_constant = sqrt(double(lpigetsol(prob,gam)));

echo off

% Print the bound and the exact value (b-a)/pi
fprintf(['\n If successful, ',num2str(poincare_constant,4),' is an upper bound on Poincare''s constant for this problem.\n'])
fprintf([' An optimal value of Poincare''s constant on s1 in [',num2str(a),',',num2str(b),'] is known to be (b-a)/pi=',num2str((b-a)/pi,4),'.\n']);
