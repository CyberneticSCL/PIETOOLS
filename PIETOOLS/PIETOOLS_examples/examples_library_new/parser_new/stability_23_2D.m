clc; clear;
%% Define dependent variables and system variable
%   PDE: x_{t} = -1/(2*lam) * x_{tt}                | lam = 1;
%                 +(C^2/2*lam)*(x_{(2,0}            | C = 1;                        Holmes, 1994 [14] Eq. (3)
%                                 + x_{(0,2)})      | ne = 1 (state size) 
%   BCs: x(s1=0) = 0,   x(s2=0) = 0,                | 
%        x(s1=1) = 0,   x(s2=1) = 0;                |
%   Use states x1 = x,    x2 = x_{t}                |
%       =>                                          |
%   PDE: x1_{t} = x2                                |
%        x2_{t} = C^2*(x1_{(2,0} + x1_{(0,2)})      |
%                   - 2*lam*x2                      |
%   BCs: x1(s1=0) = 0,   x1(s2=0) = 0,              | 
%        x1(s1=1) = 0,   x1(s2=1) = 0;              |
%        x2(s1=0) = 0,   x2(s2=0) = 0,              | 
%        x2(s1=1) = 0,   x2(s2=1) = 0;              |
pvar s t; % define independent variables
pdevar x1 x2; x1 = setInfinite(x1,2); x2 = setInfinite(x2,2);
lam=1;C = 1;
%% Define equations
%   PDE: x1_{t} = x2                                |
%        x2_{t} = C^2*(x1_{(2,0} + x1_{(0,2)}) - 2*lam*x2                      |
eq_PDE = [diff(x1,t)==x2;
                diff(x2,t)==C^2*diff(x,s1,2)+C^2*diff(x,s2,2)- 2*lam*x2]; 
%   BCs: x1(s1=0) = 0,   x1(s2=0) = 0,              | 
%        x1(s1=1) = 0,   x1(s2=1) = 0;              |
%        x2(s1=0) = 0,   x2(s2=0) = 0,              | 
%        x2(s1=1) = 0,   x2(s2=1) = 0;              |
eq_BC = [subs(x1,s1,0)==0;subs(x1,s1,1)==0;
                 subs(x1,s2,0)==0;subs(x1,s2,1)==0;
                 subs(x2,s1,0)==0;subs(x2,s1,1)==0;
                 subs(x2,s2,0)==0;subs(x2,s2,1)==0]; 
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde',[eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
