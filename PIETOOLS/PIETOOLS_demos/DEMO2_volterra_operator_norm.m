% DEMO2_volterra_operator_norm.m
% See Chapter 11.2 of the manual for a description.
%
% This document illustrates how the norm of the Volterra integral operator
% can be computed in PIETOOLS.

% This example is also included in the paper (page 5, Demonstration 2)
% link: https://arxiv.org/pdf/1910.01338.pdf

% Volterra integral operator
% T x (s) = int(x(\theta),d\theta,a,s)
% Domain:  s,theta \in [a,b] = [0,1]

% Optimization Problem (LPI)
% min γ, such that
% T^*T ≤ γ
%%
clc; clear;
echo on
%%%%%%%%%%%%%%%%%% Start Code Snippet %%%%%%%%%%%%%%%%%%

%%% Define the operator (T x)(s) = int_{a}^{s} x(r) dr on [a,b]=[0,1]
a=0;    b=1;
opvar Top;
Top.R.R1 = 1;   Top.I = [a,b];


%%% Solve the LPI Top'*Top - gam<=0
% First, define dpvar gam and set up an optimization problem
vars = [Top.var1;Top.var2];     % Free vars in optimization problem (no optimization over these vars)
dpvar gam;                      % Decision var in optimization problem (we will minimize gam)
prob = sosprogram(vars,gam);

% Next, set gam as objective function min{gam}
prob = sossetobj(prob, gam);

% Then, enforce the constraint Top'*Top-gam<=0
opts.psatz = 1;                             % Allow Top'*Top-gam>0 outside of interval [a,b]
prob = lpi_ineq(prob,-(Top'*Top-gam),opts); % lpi_ineq(prob,Q) enforces Q>=0

% Finally, solve and retrieve the solution
prob = sossolve(prob);
operator_norm = sqrt(double(sosgetsol(prob,gam)));

%%%%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
echo off

fprintf(['\n If successful, ',num2str(operator_norm),' is an upper bound on the norm of the Volterra integral operator.\n'])
fprintf([' The exact operator norm of the Volterra integral operator is 2/pi=',num2str(2/pi),'.\n']);