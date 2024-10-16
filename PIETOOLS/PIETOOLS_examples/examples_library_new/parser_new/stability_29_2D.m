clc; clear;
%% Define dependent variables and system variable
%   ODE: X_{t} = A*X(t) + A1*X(t-tau) + B*x(t,s1=0) | tau = 1;
%   PDE: x_{t} = x_{s1s1} + a*x + a2*x(t-tau);      | A = -1;    A1 = 0.5;          Kang, 2017 [19] 
%   BCs: x_{s1}(t,s1=0) = 0;                        | B = -1;
%        x(t,s1=1) = 0;                             | a = 1;    a2 = 0.5;
pvar s t; % define independent variables
odevar X; pdevar x; 
tau = 1; A = -1;    A1 = 0.5;a = 1;    a2 = 0.5;B = -1;
%% Define equations
%   ODE: X_{t} = A*X(t) + A1*X(t-tau) + B*x(t,s1=0) | tau = 1;
%   PDE: x_{t} = x_{s1s1} + a*x + a2*x(t-tau);      | A = -1;    A1 = 0.5;          Kang, 2017 [19] 
eq_ODE = diff(X,t)==A*X+A1*subs(X,t,t-tau)+B*subs(x,s,0);
eq_PDE = diff(x,t)==diff(x,s,2)+a*x+a2*subs(x,t,t-tau);
%   BCs: x_{s1}(t,s1=0) = 0;                        | B = -1;
%        x(t,s1=1) = 0;                             | a = 1;    a2 = 0.5;
eq_BC = [subs(diff(x,s),subs,0)==0;subs(x,s,1)==0];        
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde', [eq_ODE;eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
