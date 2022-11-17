clc; clear;
pvar s t theta; % define independent variables
%% Define dependent variables and system variable
% 	ODE:  xo_{t} = k*xo                             | k = -1              
%   PDE:  x_{t} = x_{ss}                            |
%   BCs:  x(s=0) = 0,     x(s=1) = xo               |
x=state('ode'); X = state('pde');
pde = sys(); 
k = -1;
%% Define equations
eq_ODE=diff(x,t)==k*x; % ODE: xo_{t} = k*xo
eq_PDE = diff(X,t)==diff(X,s,2); % 	PDE: x_{t} = x_{ss}
eq_BC = [subs(X,s,0)==0;subs(X,s,1)==x]; %  BCs: x(s=a) = 0,     x(s=b) = 0 
%% addequations to pde system; set control inputs/observed inputs, if any
pde = addequation(pde,[eq_ODE;eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
