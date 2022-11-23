clc; clear;
pvar s t theta; % define independent variables
%% Define dependent variables and system variable
%   PDE: x_{t} = x_{ss} + s*w(t)                    |                       (gamma = 0.3333)
%   BCs: x(s=0) = 0,        x_{s}(s=1) = 0          |
%   Out: z(t) = int(x(t,s),s,0,1)                   |         
x=state('pde'); z = state('out'); w=state('in');
pde = sys();     
%% Define equations
eq_PDE = [diff(x,t)==-diff(x,s,2)+s*w]; %   PDE: x_{t} = x_{ss} + s*w(t)
eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,1)==0]; %   BCs: x(s=0) = 0, x_{s}(s=1) = 0
eq_out = [z==int(x,s,[0,1])]; %   Out: z(t) = int(x(t,s),s,0,1)
%% addequations to pde system; set control inputs/observed inputs, if any
pde = addequation(pde,[eq_PDE;eq_BC;eq_out]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
