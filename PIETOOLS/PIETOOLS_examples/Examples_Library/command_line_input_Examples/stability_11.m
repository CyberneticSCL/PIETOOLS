clc; clear;
pvar s t theta; % define independent variables
%% Define dependent variables and system variable
% 	PDE:  x_{t} = Cm*x + (1/R)*x_{ss}               | R = 2.7
%   BCs:  x(s=0) = 0,     x_{s}(s=1) = 0            | Cm = [1, 1.5; 5, 0.2]     
x=state('pde',2); 
pde = sys(); 
R = 2.7;Cm = [1, 1.5; 5, 0.2];
%% Define equations
eq_PDE = diff(x,t)==Cm*x+(1/R)*diff(x,s,2); % 	PDE: x_{t} = Cm*x + (1/R)*x_{ss}
eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,1)==0]; %  BCs: x(s=0) = 0,     x_{s}(s=1) = 0 
%% addequations to pde system; set control inputs/observed inputs, if any
pde = addequation(pde,[eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
