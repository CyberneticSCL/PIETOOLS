clc; clear;
%% Define dependent variables and system variable
%   PDE: x_{t} = x_{ss} + w(t)                      | ne = 1 (state size)   (Answer: 1.0045) 
%   BCs: x(s=0) = 0,        x(s=1) = 0              |
%   Out: z(t) = int(x(t,s),s,0,1) + w(t)            |
%        y(t) = x_{s}(s=1)                          |
pvar s t; % define independent variables
pdevar x; outputvar y z; inputvar w;
%% Define equations
eq_PDE = [diff(x,t)==diff(x,s,2)+w]; %   PDE: x_{t} = x_{ss} + w(t)
eq_BC = [subs(x,s,0)==0;subs(x,s,1)==0]; %   BCs: x(s=0) = 0, x(s=1) = 0
eq_out = [z==int(x,s,[0,1])+w; y== subs(x,s,1)]; %   Out: z(t) = int(x(t,s),s,0,1) y(t) = x_{s}(s=1)
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde',[eq_PDE;eq_BC;eq_out]);
pde = setObserve(sys,y);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
