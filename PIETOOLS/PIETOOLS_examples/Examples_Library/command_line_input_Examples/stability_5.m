clc; clear; clear stateNameGenerator
pvar s t; % define independent variables
%% Define dependent variables and system variable
%   PDE: x_{t} = lam*x + x_{ss}                     | lam = 9.86
%   BCs: x(s=0) = 0,      x(s=1) = 0                |
x = pde_var(s,[0,1]);
lam = 9.86;
%% Define equations
eq_PDE = diff(x,t)==lam*x+diff(x,s,2); % 	PDE: x_{t} = lam*x + x_{ss}
eq_BC = [subs(x,s,0)==0;subs(x,s,1)==0]; %  BCs: x(s=0) = 0,      x(s=1) = 0
%% initialize pde system;
pde = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);



% %% Define dependent variables and system variable
% %   PDE: x_{t} = lam*x + x_{ss}                     | lam = 9.86
% %   BCs: x(s=0) = 0,      x(s=1) = 0                |
% x = state('pde');
% pde = sys(); lam = 9.86;
% %% Define equations
% eq_PDE = diff(x,t)==lam*x+diff(x,s,2); % 	PDE: x_{t} = lam*x + x_{ss}
% eq_BC = [subs(x,s,0)==0;subs(x,s,1)==0]; %  BCs: x(s=0) = 0,      x(s=1) = 0
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
