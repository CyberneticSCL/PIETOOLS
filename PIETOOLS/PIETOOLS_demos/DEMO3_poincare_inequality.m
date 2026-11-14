% DEMO3_poincare_inequality.m
% See Chapter 11.3 of the manual for a description.
%
% This document illustrates how the Poincare constant can be found
% using PIETOOLS.
% The Pointcare constant is the smallest value c such that
%   ||x|| ≤ c||x_s||    for all x in {x: x_{s}\in L_2[0,1] & x(0)=x(1)=0}
% We can pose this as an optimization program
%   min_{c}     c^2,
%      s.t.     <x, x> − c^2 <x_s,x_s> ≤ 0
%               x \in H1:={x: x_{s}\in L_2[0,1] & x(0)=x(1)=0}
% To solve, we suppose that x is twice-differentiable,
%               D := {x: x_{ss}\in L_2[0,1] & x(0)=x(1)=0}.
% We define a boundary value problem
%               x_{ss} = v,     x(0)=x(1)=0,    v in L2[0,1]
% and find the solution H2op: L2 --> D such that
%               x = H2op*v
% for v in L2[0,1] free of boundary conditions. We also compute 
% H1op = (d/ds) H2op such that
%     (H1op*v)(s) = x_{s}(s).
% Given these operators, we can pose the Poincare optimization problem as
% an LPI
%   min_{gam}   gam,
%      s.t.     H2op'*H2op -gam*H1op'*H1op ≤ 0,
% where now c = sqrt(gam).
%
% This example is also included in the paper (page 6, Demonstration 3)
% link: https://arxiv.org/pdf/1910.01338.pdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - DEMO3
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
% MP, SS, DJ, 2022: Initial coding;
% DJ, 10/20/2024: Update to use new LPI programming functions;
% DJ, 12/23/2024: Simplify demo (remove lines of code where possible);
% DJ (CC), 10/30/2026: Declare the LPI using the sopvar methods;
% MMP, 11/04/2026: The LPI program is declared by LPIPROGRAM on the spatial
%   variables of the PIE; LPIPROGRAM_SOP, which read them off T, is retired.
% DJ, 11/06/2026: Use BVP instead of PDE
% DJ, 11/11/2026: Use `convert` to solve the BVP;

clc; clear; clear stateNameGenerator;
echo on

% =============================================
% === Declare the operators of interest

% % Declare the BVP
% %     x_{ss} = v,     x(a) = x(b) = 0
a = 0;  b = 1;
pvar t s
x = pde_var(s,[a,b],2);     % x in W^2[a,b]
BVP = [subs(x,s,a)==0;
       subs(x,s,b)==0];

% % Solve the boundary value problem
% %     x = Top * v     v in L2[0,1]
BVP_sol = convert(BVP);
H2op = BVP_sol.T;                   % (H2op*x_{ss}) = x;
H1op = diff(H2op,s);                % (H1op*x_{ss}) = x_{s}


% =============================================
% === Declare the LPI

% % Initialize LPI program
prob = lpiprogram(H2op.vars,H2op.dom);

% % Declare decision variables:
% %   gam \in \R
[prob,gam] = lpidecvar(prob,'gam');     % scalar decision variable

% % Set inequality constraints:
% %   gam*H1op'*H1op' - H2op'*H2op >= 0
% % (by default, this is allowed to fail outside of [a,b])
prob = lpi_ineq_sop(prob,gam*(H1op'*H1op)-H2op'*H2op);                      % DJ (CC), 10/30/2026

% % Set objective function:
% %   min gam
prob = lpisetobj(prob, gam);

% % Solve LPI and retrieve solution
prob = lpisolve(prob);
poincare_constant = sqrt(double(lpigetsol_sop(prob,gam)));                  % DJ (CC), 10/30/2026


echo off

fprintf(['\n If successful, ',num2str(poincare_constant,4),' is an upper bound on Poincare''s constant for this problem.\n'])
fprintf([' An optimal value of Poincare''s constant on domain [0,1] is known to be 1/pi=',num2str(1/(pi),4),'.\n']);