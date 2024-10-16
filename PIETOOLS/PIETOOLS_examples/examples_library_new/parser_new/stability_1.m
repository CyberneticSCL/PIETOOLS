clc; clear;
%% Define dependent variables and system variable
% 	PDE: x_{t} = v*x_{s}                                 
%   BCs: x(s=0) = 0 
pvar s t; % define independent variables
pdevar x; v=-1;
%% Define equations
eq_PDE = diff(x,t)==v*diff(x,s); % 	PDE: x_{t} = x_{s}
eq_BC = subs(x,s,0)==0;        %   BCs: x(s=0) = 0 
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde', [eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
