clc; clear;
%% Define dependent variables and system variable
%   PDE: xi_{t} = w(t) + lamb*xi        i=1:ne      | lam = (1-1e-2)*pi^2   (gamma = 8.1069 for lamb = (1-1e-2)*pi^2)         
%                 + sum(xk_{ss},k=1,i)              | ne = 1 (state size)           Shivakumar 2019 [12] (Example 3)
%   BCs: x(s=0) = 0,        x(s=1) = 0              |   
%   Out: z(t) = int(x(t,s),s,0,1)                 |         
pvar s t; % define independent variables
pdevar x; outputvar z; inputvar w;
ne = 1; x.len = ne; z.len = ne;
%% Define equations
%   PDE: xi_{t} = w(t) + lamb*xi        i=1:ne      | lam = (1-1e-2)*pi^2   (gamma = 8.1069 for lamb = (1-1e-2)*pi^2)         
%                 + sum(xk_{ss},k=1,i)              | ne = 1 (state size)           Shivakumar 2019 [12] (Example 3)
I = triu(ones(ne))+eye(ne);
eq_PDE = [diff(x,t)==I*diff(x,s,2)+lamb*x+w]; 
eq_BC = [subs(x,s,0)==0;subs(x,s,1)==0]; %   BCs: x(s=0) = 0, x_{s}(s=1) = 0
eq_out = [z==int(x,s,[0,1])]; %   Out: z(t) = int(x(t,s),s,0,1)
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde',[eq_PDE;eq_BC;eq_out]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
