clc; clear; clear stateNameGenerator
pvar s t; % define independent variables
%     | PDE: u_{tt} + 2*ad*u_{t} = -ad^2*u + u_{ss} | k = 1                         Datko 1986 [9] (Test 7.5d)
%     | BCs: u(s=0) = 0                             |
%     |      u_{s}(s=1) = -k*u_{t}(s=1)             |
%
% use states x1 = u_{t}, x2 = u, x3 = u_{s}.
% % Then 
% x1_{t} = -2*ad*x1 - ad^2*x2 + x3_{s}       |
% x2_{t} = x1 
% % x3_{t} = x1_{s} 
% x1(0) = 0, x2(0) = 0,  k*x1(1) + x3(1) = 0.
k = 1;    ad = 1;
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
pde_var x1 x2 x3
x1.vars = s;    x2.vars = s;  x3.vars = s; 
%% Define equations
eq_PDE = [diff(x1,t,1)==-2*ad*x1-ad^2*x2+diff(x3,s,1); %   PDE: x1_{t} = -2*ad*x1 - ad^2*x2 + x3_{s}       |
                  diff(x2,t,1)==x1;% x2_{t} = x1 
                  diff(x3,t,1)==diff(x1,s,1)];% x3_{t} = x1_{s} 
eq_BC = [subs(x1,s,0)==0;subs(x2,s,0)==0;%%   BCs: x1(0) = 0,     x2(0) = 0,                  |                  |
                k*subs(x1,s,1)+subs(x3,s,1)==0;]; %  k*x1(1) + x3(1) = 0     
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);