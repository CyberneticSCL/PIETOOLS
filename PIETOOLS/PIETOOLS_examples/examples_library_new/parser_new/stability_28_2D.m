clc; clear;
%% Define dependent variables and system variable
%   PDE: x_{t} = x_{s1s1} + lam * x;    s1 in [0,1] | tau = 1;
%   BCs: x(t,s1=0) = 0;                             | lam = 1;                      Kristic, 2009 [18] 
%        x(t,s1=1) = u(t-tau);                      |
pvar s t; % define independent variables
pdevar x; inputvar u;
tau = 1; lam = 1;
%% Define equations
eq_PDE = diff(x,t)==diff(x,s,2)+lam*x ; %   PDE: x_{t} = x_{s1s1} + lam * x;
%   BCs: x(t,s1=0) = 0;                             | lam = 1;                      Kristic, 2009 [18] 
%        x(t,s1=1) = u(t-tau);                      |
eq_BC = [subs(x,s,0)==0;subs(x,s,1)==subs(u,t,t-tau)];        
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde', [eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
