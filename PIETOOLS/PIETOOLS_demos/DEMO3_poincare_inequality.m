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
%               {x: x_{ss}\in L_2[0,1] & x(0)=x(1)=0}.
% We define an artificial PDE on x as
%       d/dt x(t,s) = x_{ss}(t,s)
%            z(t,s) = x_{s}(t,s)
%            x(t,0) = x(t,1) = 0
% For this PDE, we can obtain a PIE representation
%       d/dt (H2op*v)(t,s) = v(t,s)
%                   z(t,s) = (H1op*v)(t,s)
% where now v(t) is free of boundary conditions, and H1op and H2op are s.t.
%   (H2op*x_{ss})(s) = x(s),
%   (H1op*x_{ss})(s) = x_{s}(s).
% Given these operators, we can pose the Poincare optimization problem as
% an LPI
%   min_{gam}   gam,
%      s.t.     H2op'*H2op -gam*H1op'*H1op ≤ 0,
% where now c = sqrt(gam).
%
% This example is also included in the paper (page 6, Demoenstration 3)
% link: https://arxiv.org/pdf/1910.01338.pdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - DEMO3
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
% MP, SS, DJ, 2022: Initial coding;
% DJ, 10/20/2024: Update to use new LPI programming functions;
% DJ, 12/23/2024: Simplify demo (remove lines of code where possible);

clc; clear; clear stateNameGenerator;
echo on

% =============================================
% === Declare the operators of interest

% % Declare system as PDE
a = 0;  b = 1;
pvar t s
x = pde_var(1,s,[a,b]);
z = pde_var('output',1,s,[a,b]);
PDE = [diff(x,t)==diff(x,s,2);
       z==diff(x,s,1);
       subs(x,s,a)==0;
       subs(x,s,b)==0];

% % Convert PDE to PIE
PIE = convert(PDE);
H2op = PIE.T;       % (H2op*x_{ss}) = x;
H1op = PIE.C1;      % (H1op*x_{ss}) = x_{s}


% =============================================
% === Declare the LPI

% % Initialize LPI program
prob = lpiprogram(s,[a,b]);

% % Declare decision variables:
% %   gam \in \R
[prob,gam] = lpidecvar(prob,'gam');     % scalar decision variable

% % Set inequality constraints:
% %   gam*H1op'*H1op' - H2op'*H2op >= 0
opts.psatz = 1;                 % allow gam H1op'*H1op < H2op'*H2op outside of [a,b]
prob = lpi_ineq(prob,gam*(H1op'*H1op)-H2op'*H2op,opts);

% % Set objective function:
% %   min gam
prob = lpisetobj(prob, gam);

% % Solve LPI and retrieve solution
prob = lpisolve(prob);
poincare_constant = sqrt(double(lpigetsol(prob,gam)));


echo off

fprintf(['\n If successful, ',num2str(poincare_constant,4),' is an upper bound on Poincare''s constant for this problem.\n'])
fprintf([' An optimal value of Poincare''s constant on domain [0,1] is known to be 1/pi=',num2str(1/(pi),4),'.\n']);