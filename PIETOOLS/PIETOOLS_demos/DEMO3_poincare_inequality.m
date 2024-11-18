% DEMO3_poincare_inequality.m
% See Chapter 11.3 of the manual for a description.
%
% This document illustrates how the Poincare constant can be found
% using PIETOOLS.
% The Pointcare constant is the smallest value c such that
%   ||x|| ≤ c||x_s||    for all x in {x: x_{s}\in L_2[0,1] & x(0)=x(1)=0}
% We can pose this as an optimization program
%   min   c,
%   s.t.  <x, x> − c <x_s,x_s> ≤ 0
%           x in {x: x_{s}\in L_2[0,1] & x(0)=x(1)=0}
% To solve, we use PIETOOLS PDE to PIE converter to compute H1op and H2op
% such that for all x satisfying x(0)=x(1)=0,
%   x = H2op*x_{ss}     and     x_{s} = H1op*x_{s};
% Then, we pose the Poincare inequality as an LPI
%   min   c,
%   s.t.  H2'*H2 -c*H1'*H1 ≤ 0

% This example is also included in the paper (page 6, Demoenstration 3)
% link: https://arxiv.org/pdf/1910.01338.pdf

clc; clear; clear stateNameGenerator;
echo on
%% %%%%%%%%%%%%%%%% Start Code Snippet %%%%%%%%%%%%%%%%%%

% =============================================
% === Declare the operators of interest

% Declare a PDE 
%   d/dt x = x_{ss}
%       z1 = x_{s}
%   x(s=a) = x(s=b) = 0;
pvar s;
a = 0;      b = 1;
x = pde_var('state',1,s,[a,b]);
z = pde_var('out',1,s,[a,b]);
PDE = [diff(x,'t')==diff(x,s,2);
       z==diff(x,s);
       subs(x,s,a)==0;
       subs(x,s,b)==0];

% Convert to a PIE and extract operators
%   d/dt H2op*x_{ss} = x_{ss};
%                 z1 = H1op*x_{ss};
PIE = convert(PDE);
H2op = PIE.T;       % H2op*x_{ss} = x
H1op = PIE.C1;      % H1op*x_{ss} = x_{s}

% =============================================
% === Declare the LPI

% Initialize LPI program
prob = lpiprogram(PIE.vars,PIE.dom);

% Set inequality constraints:
%   H2'*H2-gam H1'*H1 <= 0
dpvar gam
prob = lpidecvar(prob,gam);
opts.psatz = 1;                 % Allow H2'*H2 > gam H1'*H1 outside of [a,b]
prob = lpi_ineq(prob,-(H2op'*H2op-gam*H1op'*H1op),opts);

% Set objective function:
%   min gam
prob = lpisetobj(prob, gam);

% Solve and retrieve the solution
prob = lpisolve(prob);
poincare_constant = sqrt(double(sosgetsol(prob,gam)));

%% %%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
echo off

fprintf(['\n If successful, ',num2str(poincare_constant),' is an upper bound on Poincare''s constant for this problem.\n'])
fprintf([' An optimal value of Poincare''s constant on domain [0,1] is known to be 1/pi=',num2str(1/(pi)),'.\n']);