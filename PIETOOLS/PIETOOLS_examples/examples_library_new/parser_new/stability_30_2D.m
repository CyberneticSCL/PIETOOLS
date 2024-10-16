clc; clear;
%% Define dependent variables and system variable
%   PDE:    x_{t}  = c*x_{s1s1}(t,s1) + a0*x(t,s1)  | tau = 1;
%                       - a1*x(t-tau,s1);           | c = 1;
%   BCs:    x(t,s1=0) = 0;                          | a0 = 1.9;
%           x(t,s1=pi) = 0;                         | a1 = 1;
pvar s t; % define independent variables
 pdevar x; x.dom = [0,pi];
tau = 1; c = 1;a0 = 1.9;a1 = 1;
%% Define equations
%   PDE:    x_{t}  = c*x_{s1s1}(t,s1) + a0*x(t,s1)  | tau = 1;
%                       - a1*x(t-tau,s1);           | c = 1;
eq_PDE = diff(x,t)==c*diff(x,s,2)+a0*x-a1*subs(x,t,t-tau);
%   BCs:    x(t,s1=0) = 0;                          | a0 = 1.9;
%           x(t,s1=pi) = 0;                         | a1 = 1;
eq_BC = [subs(x,s,0)==0;subs(x,s,pi)==0];        
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde', [eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
