clc; clear;
%% Define dependent variables and system variable
%   PDE: phi_{tt} = phi_{ss} + w(t)                 | k = 0.5               (gamma = 2 for k = 0.5)  
%   BCs: phi(s=0) = 0                               |
%        phi_{s}(s=1) = -k*phi_{t}(s=1)             |
%   Out: z(t) = int(phi_{t}(t,s),s,0,1)             |
%   Use states x1 = phi_{s},    x2 = phi_{t}        |
%       =>                                          |
%   PDE: x1_{t} = x2_{s}                            |
%        x2_{t} = x1_{s} + w(t)                     |
%   BCs: x2(0) = 0,     x1(1) + k*x2(1) = 0         |
pvar s t; % define independent variables
pdevar x1 x2; outputvar z; inputvar w;
k = 0.5
%% Define equations
%   PDE: x1_{t} = x2_{s}                            |
%        x2_{t} = x1_{s} + w(t)                     |
eq_PDE = [diff(x1,t)==diff(x2,s);
                  diff(x2,t)==diff(x1,s)+w];
%   BCs: x2(0) = 0,     x1(1) + k*x2(1) = 0         |
eq_BC = [subs(x2,s,0)==0; subs(x1,s,1) + k*subs(x2,s,1) = 0]; 
eq_out = [z==int(x2,s,[0,1])]; %   Out: z(t) = int(x2(t,s),s,0,1)
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde',[eq_PDE;eq_BC;eq_out]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
