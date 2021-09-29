% This document has an illustration of finding the poincare constant

% This example is included in the paper (page 6, Demoenstration 3)
% link: https://arxiv.org/pdf/1910.01338.pdf


% What is Poincare Inequality?
% Find C; such that for an function u \in H^1[0, 1]
% ||u|| ≤ C||u_s||

% Optimization Problem
% min C, such that
% <u, u> − C <u_s,u_s> ≤ 0

%%
%%%%%%%%%%%%%%%%%% Code Snippet %%%%%%%%%%%%%%%%%%


%%% Zero step: Define the polynomial variable s, if needed
clc; clear; 
pvar s theta;

%%% First step: Define the Domain [a,b]
PDE.dom = [0,1];

%%% Second step: u is twice differentiated and hence is in n2
% number of pde states are defined as n2=1 and n0=n1=0
PDE.n0=0; PDE.n1=0; PDE.n2=1;

%%% Make sure the first derivative of u is computed as: u_s = PIE.A u_ss
PDE.A1 = 1;

%%% Define Boundary conditions in the form B[xp(a); xp(b); xp_s(a) xp_s(b)]= 0
PDE.B = [1 0 0 0; 0 1 0 0]; 

%%% Solve the LPI < Top u_ss, Top u_ss> - gam <H1 u_ss, H1 u_ss> <=0
%%% where (Top u_ss) = u and (H1 u_ss) = u_s
% First, define pvar gam and set up an optimization problem
pvar gam;
prob = sosprogram([s;theta],gam);

%%% Use conversion script to convert 
% This defines the PI operators Top and H1
PIE = convert_PIETOOLS_PDE_batch(PDE);
Top = PIE.T;
H1 = PIE.A;

% Set gam as objective function
prob = sossetobj(prob, gam);

% Since we need the negative definiteness on the interval [0,1], lets add
% psatz term
opts.psatz = 1;
% Set up the constraint Top'*Top-gam<=0
prob = lpi_ineq(prob,-(Top'*Top-gam*H1'*H1),opts);

% Solve and retrieve the solution
prob = sossolve(prob);
poincare_constant = sqrt(double(sosgetsol(prob,gam)));