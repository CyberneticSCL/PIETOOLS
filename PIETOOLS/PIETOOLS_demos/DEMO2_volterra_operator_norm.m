% DEMO2_volterra_operator_norm.m
% See Chapter 11.2 of the manual for a description.
%
% This document illustrates how an upper bound on the norm of the Volterra 
% integral operator can be computed in PIETOOLS.
% Volterra integral operator
%   (Top*x)(s) = int(x(r),dr,a,s),        s in [a,b]
% Optimization Problem (LPI)
% min   γ,
% s.t   T'*T ≤ γ

% This example is also included in the paper (page 5, Demonstration 2)
% link: https://arxiv.org/pdf/1910.01338.pdf

clc; clear; clear stateNameGenerator;
echo on

% =============================================
% === Declare the operators of interest

% % Define Volterra operator
% %  (Top*x)(s) = int_{a}^{s} x(r) dr    s in [a,b]
opvar Top;
a=0;    b=1;    Top.I = [a,b];
Top.R.R1 = 1;                   

% =============================================
% === Declare the LPI

% % Initialize LPI program
prob = lpiprogram(Top.vars,Top.I);

% % Set inequality constraints
% %   Top'*Top-gam <= 0
dpvar gam
prob = lpidecvar(prob,gam);
opts.psatz = 1;                             % Allow Top'*Top-gam>0 outside of [a,b]
prob = lpi_ineq(prob,-(Top'*Top-gam),opts); % lpi_ineq(prob,Q) enforces Q>=0

% % Set objective function:
% %   min gam
prob = lpisetobj(prob, gam);

% % Solve LPI and retrieve solution
prob = lpisolve(prob);
operator_norm = sqrt(double(lpigetsol(prob,gam)));


echo off

fprintf(['\n If successful, ',num2str(operator_norm),' is an upper bound on the norm of the Volterra integral operator.\n'])
fprintf([' The exact operator norm of the Volterra integral operator is 2/pi=',num2str(2/pi),'.\n']);