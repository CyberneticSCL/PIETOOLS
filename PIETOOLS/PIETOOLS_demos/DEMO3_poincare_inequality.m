% DEMO3_poincare_inequality.m
% See Chapter 11.3 of the manual for a description.
%
% This document illustrates how the Poincare constant can be found
% using PIETOOLS.
% The Pointcare constant is the smallest value c such that
%   ||x|| ≤ c||x_s||    for all x in {x: x_{s}\in L_2[0,1] & x(0)=x(1)=0}
% We can pose this as an optimization program
%   min_{c}     c,
%      s.t.     <x, x> − c <x_s,x_s> ≤ 0
%               x \in H1:={x: x_{s}\in L_2[0,1] & x(0)=x(1)=0}
% To solve, define operators
%   (Hop*v)(s) := int_{0}^{s} v(r) dr,
%    Iop*v     := int_{0}^{1} v(s) ds.
% Then, for all x \in H1, we have x = Hop*x_{s} and Iop*x_{s}=0.
% Conversely, for all v satisfying Iop*v=0, we have Hop*v \in H1.
% Then, using Finsler's lemma, we can pose the problem as an LPI
%   min_{c,Xop} c,
%      s.t.     Hop'*Hop -c +Xop'*Iop +Iop*Xop' ≤ 0.

% This example is also included in the paper (page 6, Demoenstration 3)
% link: https://arxiv.org/pdf/1910.01338.pdf

clc; clear; clear stateNameGenerator;
echo on
%% %%%%%%%%%%%%%%%% Start Code Snippet %%%%%%%%%%%%%%%%%%

% =============================================
% === Declare the operators of interest

% Declare Hop:L2-->L2 and Iop:L2-->R
opvar Hop Iop;
a = 0;  b = 1;     
Hop.I = [a,b];  Iop.I = [a,b];
Hop.R.R1 = 1;       % (Hop*v)(s) = int_{a}^{s} v(r)dr
Iop.Q1 = 1;         % (Iop*v) = int_{a}^{b} v(s)ds


% =============================================
% === Declare the LPI

% Initialize LPI program
prob = lpiprogram(Hop.vars,Hop.dom);

% Declare decision variables:
%   gam \in \R   and    Xop:L2-->\R
dpvar gam                               % scalar decision variable
prob = lpidecvar(prob,gam);
[prob,Xop] = lpivar(prob,Iop.dim,5);    % operator decision variable

% Set inequality constraints:
%   Hop'*Hop -gam  +Xop'*Iop +Iop*Xop' <= 0
opts.psatz = 1;                 % Allow Hop'*Hop > gam outside of [a,b]
prob = lpi_ineq(prob,-(Hop'*Hop -gam +Xop'*Iop +Iop'*Xop),opts);

% Set objective function:
%   min gam
prob = lpisetobj(prob, gam);

% Solve and retrieve the solution
prob = lpisolve(prob);
poincare_constant = sqrt(double(lpigetsol(prob,gam)));


%% %%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
echo off

fprintf(['\n If successful, ',num2str(poincare_constant,7),' is an upper bound on Poincare''s constant for this problem.\n'])
fprintf([' An optimal value of Poincare''s constant on domain [0,1] is known to be 1/pi=',num2str(1/(pi),7),'.\n']);