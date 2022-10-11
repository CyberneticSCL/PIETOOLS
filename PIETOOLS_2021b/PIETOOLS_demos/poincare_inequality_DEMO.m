% This document has illustrates how the poincare constant can be found
% using PIETOOLS

% This example is included in the paper (page 6, Demoenstration 3)
% link: https://arxiv.org/pdf/1910.01338.pdf


% What is Poincare Inequality?
% Find C; such that for an function u \in H^1[0, 1]
% ||u|| ≤ C||u_s||

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
Top = PIE.T;
H1 = PIE.A;


%%%%%   Solve the LPI < Top u_ss, Top u_ss> - gam <H1 u_ss, H1 u_ss> <=0
%%%%%   where (Top u_ss) = u and (H1 u_ss) = u_s
% % First, define dpvar gam and set up an optimization problem
dpvar gam;
vars = [Top.var1; Top.var2];
prob = sosprogram(vars,gam);

% % Set gam as objective function
prob = sossetobj(prob, gam);

% % Set up the constraint Top'*Top-gam<=0
opts.psatz = 1;     % Add psatz term to allow Top'*Top>gam outside of [a,b]
prob = lpi_ineq(prob,-(Top'*Top-gam*H1'*H1),opts);

% Solve and retrieve the solution
prob = sossolve(prob);
poincare_constant = sqrt(double(sosgetsol(prob,gam)));

%%%%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
echo off

fprintf(['\n If successful, ',num2str(poincare_constant),' is an upper bound on Poincare''s constant for this problem.\n'])
fprintf([' An optimal value of Poincare''s constant for this problem is known to be 1/pi=',num2str(1/(pi)),'.\n']);