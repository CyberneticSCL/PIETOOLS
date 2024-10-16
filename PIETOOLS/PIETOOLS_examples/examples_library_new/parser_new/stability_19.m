clc; clear;
%% Define dependent variables and system variable
%   Use states x1 = u_{t}, x2 = u, x3 = u_{s}.      | k = 1                
%       =>                                          | ad = 1
%   PDE: x1_{t} = -2*ad*x1 - ad^2*x2 + x3_{s}       |
%   BCs: x1(0) = 0,     x2(0) = 0,                  |
%        k*x1(1) + x3(1) = 0                        |
pvar s t; % define independent variables
pdevar x1 x2 x3;
k=1;ad = 1;
%% Define equations
%   PDE: x1_{t} = -2*ad*x1 - ad^2*x2 + x3_{s}       |
%            x2_{t} = x1    |
%            x3_{t} = x2_{s}
eq_PDE = [diff(x1,t)==-2*ad*x1 - ad^2*x2 + diff(x3,s,2);
                  diff(x2,t)==x1;diff(x3,t)==diff(x2,s)];
%   BCs: x1(0) = 0,     x2(0) = 0,                  |
%        k*x1(1) + x3(1) = 0                        |
eq_BC = [subs(x1,s,0)==0; subs(x2,s,0)==0;
               k*subs(x1,s,1)+subs(x3,s,1)==0]; 
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde',[eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
