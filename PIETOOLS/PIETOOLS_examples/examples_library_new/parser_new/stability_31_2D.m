clc; clear;
%% Define dependent variables and system variable
%   PDE:    x_{tt}  = x_{ss}(t,s);   s in [0,1]     | tau = 1;  
%   BCs:    x(t,0) = 0;                             | k = 1;
%           x_{s}(t,1) = -k*(1-mu)*x_{t}(t,1)       | mu = 0.4;
%                           - k*mu*x_{t}(t-tau,1);  |
pvar s t; % define independent variables
pdevar x; 
tau = 1;k = 1;mu = 0.4;
%% Define equations
%   PDE:    x_{tt}  = x_{ss}(t,s);   s in [0,1]     | tau = 1;  
eq_PDE = diff(x,t,2)==diff(x,s,2);
%   BCs:    x(t,0) = 0;                             | k = 1;
%           x_{s}(t,1) = -k*(1-mu)*x_{t}(t,1)       | mu = 0.4;
%                           - k*mu*x_{t}(t-tau,1);  |
eq_BC = [subs(x,s,0)==0;subs(x,s,1)==-k*(1-mu)*subs(diff(x,t),s,1)-k*mu*subs(subs(diff(x,t),t,t-tau),s,1)];        
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde', [eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
