clc; clear;
%% Define dependent variables and system variable
%   ODE: xo_{t} = xo + x_{s}(s=0)               	| k = -2   
%   PDE: x_{t} = x_{ss}                             |                            
%   BCs: x(s=0) = -xo,    x(s=1) = k*xo             |
pvar s t; % define independent variables
odevar x; pdevar X;
k = -2;
%% Define equations
eq_ODE=diff(x,t)==x+subs(diff(X,s),s,0); % ODE: xo_{t} = xo + x_{s}(s=0)
eq_PDE = diff(X,t)==diff(X,s,2); % 	PDE: x_{t} = x_{ss}
eq_BC = [subs(X,s,0)==-x;subs(X,s,1)==k*x]; %  BCs: x(s=0) = -xo,    x(s=1) = k*xo
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde',[eq_ODE;eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
