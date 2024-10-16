clc; clear;
%% Define dependent variables and system variable
%   PDE: u_{tt} = u_{ss}                            | k=1                   (stable for k=1)                  
%   BCs: u(s=0) = 0,                                |                               Peet 2019 [8] (Example 8.2)
%        u_{s}(s=1) = -k*u_{t}(s=1)                 |
%   Use states x1 = u_{s}, x2 = u_{t}.              |
%       =>                                          |
%   PDE: x1_{t} = x2_{s},   x2_{t} = x1_{s}         |
%   BCs: x2(0) = 0,         x1(1) + k*x2(1) = 0     |
pvar s t; % define independent variables
pdevar x1 x2;
k=1;
%% Define equations
%   PDE: x1_{t} = x2_{s},   x2_{t} = x1_{s}         |
eq_PDE = [diff(x1,t)==diff(x2,s); 
          diff(x2,t)==diff(x1,s)];
%   BCs: x2(0) = 0,         x1(1) + k*x2(1) = 0     |
eq_BC = [subs(x2,s,0)==0;subs(x1,s,1)+k*subs(x2,s,1)==0]; 
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde',[eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
