% This document has an illustration of finding the operator norm of
% volterra integral operator

% This example is included in the paper (page 5, Demoenstration 2)
% link: https://arxiv.org/pdf/1910.01338.pdf

% Volterra integral operator
% T x (s) = int(x(\theta),d\theta,a,s)
% Domain:  a=0, b=1

% Optimization Problem
% min γ, such that
% T^*T ≤ γ
%%
%%%%%%%%%%%%%%%%%% Code Snippet %%%%%%%%%%%%%%%%%%

%%% Zero step: Define the polynomial variable s, if needed
clc; clear;
pvar s theta;

%%% First step: Define the Domain [a,b]
a=0;b=1;

%%% Define the operator
opvar Top;
Top.R.R1 = 1; Top.I = [a,b];

%%% Solve the LPI Top'*Top - gam<=0
% First, define pvar gam and set up an optimization problem
pvar gam;
prob = sosprogram([s;theta],gam);

% Set gam as objective function
prob = sossetobj(prob, gam);
opts.psatz = 1;

% Set up the constraint Top'*Top-gam<=0
prob = lpi_ineq(prob,-(Top'*Top-gam),opts);

% Solve and retrieve the solution
prob = sossolve(prob);
operator_norm = sqrt(double(sosgetsol(prob,gam)));