clc; clear; clear stateNameGenerator
pvar s t; % define independent variables
%   PDE: u_{tt} = u_{ss}                            | k=1                   (stable for k=1)                  
%   BCs: u(s=0) = 0,                                |                               Peet 2019 [8] (Example 8.2)
%        u_{s}(s=1) = -k*u_{t}(s=1)                 |
%   Use states x1 = u_{s}, x2 = u_{t}.              |
%       =>                                          |
%   PDE: x1_{t} = x2_{s},   x2_{t} = x1_{s}         |
%   BCs: x2(0) = 0,         x1(1) + k*x2(1) = 0     |
k = 1;
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
pde_var x1 x2
x1.vars = s;    x2.vars = s; 
%% Define equations
eq_PDE = [diff(x1,t,1)==diff(x2,s,1); 
                  diff(x2,t,1)==diff(x1,s,1);];
eq_BC = [subs(x2,s,0)==0;subs(x1,s,1)+k*subs(x2,s,1)==0;] 
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);