% DEMO3_poincare_inequality.m
% See Chapter 11.3 of the manual for a description.
%
% This document illustrates how the Poincare constant can be found
% using PIETOOLS

% This example is also included in the paper (page 6, Demoenstration 3)
% link: https://arxiv.org/pdf/1910.01338.pdf


% What is Poincare Inequality?
% Find C; such that for an function u \in H^1[0, 1]
% ||u|| ≤ C||u_s||
% where H^1[0,1] := {u: u_{s}\in L_2[0,1] & u(0)=u(1)=0}

% Optimization Problem
% min C, such that
% <u, u> − C <u_s,u_s> ≤ 0

%%
clc; clear; 
echo on
%%%%%%%%%%%%%%%%%% Start Code Snippet %%%%%%%%%%%%%%%%%%

%%%%%   Declare the PDE

% % Initialize the PDE structure and spatial variable s in [a,b]
pvar s theta;
pde_struct PDE;
a = 0;      b = 1;

% % Declare the state variables x(t,s)
PDE.x{1}.vars = s;
PDE.x{1}.dom = [a,b];
PDE.x{1}.diff = 2;          % Let x be second order differentiable wrt s.

% % Declare the PDE \dot{x}(t,s) = \partial_{s} x(t,s)
PDE.x{1}.term{1}.D = 1;     % Order of the derivative wrt s

% % Declare the boundary conditions x(t,a) = x(t,b) = 0
PDE.BC{1}.term{1}.loc = a;      % Evaluate x at s=a
PDE.BC{2}.term{1}.loc = b;      % Evaluate x at s=b

% % Initialize the system
PDE = initialize(PDE);


%%%%% Convert the PDE to a PIE
PIE = convert(PDE,'pie');
H2 = PIE.T;     % H2 x_{ss} = x
H1 = PIE.A;     % H1 x_{ss} = x_{s}


%%%%%   Solve the LPI < H2 x_ss, H2 x_ss> - gam <H1 x_ss, H1 x_ss> <=0
%%%%%   where (H2 x_ss) = x and (H1 x_ss) = x_s
% % First, define dpvar gam and set up an optimization problem
dpvar gam;
vars = [H2.var1; H2.var2];
prob = sosprogram(vars,gam);

% % Set gam as objective function to minimize
prob = sossetobj(prob, gam);

% % Set up the constraint H2'*H2-gam H1'*H1<=0
opts.psatz = 1;     % Add psatz term to allow H2'*H2 > gam H1'*H1 outside of [a,b]
prob = lpi_ineq(prob,-(H2'*H2-gam*H1'*H1),opts);

% Solve and retrieve the solution
prob = sossolve(prob);
poincare_constant = sqrt(double(sosgetsol(prob,gam)));

%%%%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
echo off

fprintf(['\n If successful, ',num2str(poincare_constant),' is an upper bound on Poincare''s constant for this problem.\n'])
fprintf([' An optimal value of Poincare''s constant on domain [0,1] is known to be 1/pi=',num2str(1/(pi)),'.\n']);