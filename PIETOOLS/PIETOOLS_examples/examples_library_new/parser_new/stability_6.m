clc; clear;

%% Define dependent variables and system variable
%   | PDE: x_{t} = lam*x + x_{ss}                   | lam = 2.466 
%   | BCs: x(s=0) = 0,      x_{s}(s=1) = 0
pvar s t; % define independent variables
pdevar x;
lam = 2.466;
%% Define equations
eq_PDE = diff(x,t)==lam*x+diff(x,s,2); % 	PDE: x_{t} = lam*x + x_{ss}
eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,1)==0]; %  BCs: x(s=0) = 0,  x_{s}(s=1) = 0
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde',[eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
