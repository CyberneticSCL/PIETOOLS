clc; clear;
%% Define dependent variables and system variable
%     | PDE: u_{tt} + 2*ad*u_{t} = -ad^2*u + u_{ss} | k = 1                         Datko 1986 [9] (Test 7.5c)
%     | BCs: u(s=0) = 0                             |
%     |      u_{s}(s=1) = -k*u_{t}(s=1)             |
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
%   Use states x1 = u_{t}, x2 = u.                  | k = 1
%       =>                                          | ad = 1
%   PDE: x1_{t} = -2*ad*x1 - ad^2*x2 + x2_{ss}      |
%   BCs: x1(0) = 0,     x2(0) = 0,                  |
%        k*x1(1) + x2_{s}(1) = 0                    |
pvar s t; % define independent variables
pdevar x1 x2;
k=1;ad = 1;
%% Define equations
%   PDE: x1_{t} = -2*ad*x1 - ad^2*x2 + x2_{ss}      |
%            x2_{t} = x1    |
eq_PDE = [diff(x1,t)==-2*ad*x1 - ad^2*x2 + diff(x2,s,2);
          diff(x2,t)==x1];
%   BCs: x1(0) = 0,     x2(0) = 0,                  |
%        k*x1(1) + x2_{s}(1) = 0                    |
eq_BC = [subs(x1,s,0)==0; subs(x2,s,0)==0;
               k*subs(x1,s,1)+subs(x2,s,1)==0]; 
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde',[eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
