clc; clear;
%% Define dependent variables and system variable
%   PDE: x_{tt} = c1*x_{(2,0)} + c2*x_{(0,2)}       | c1 = 1;   c2 = 1;
%   BCs: x(s1=0) = 0;       x(s2=0) = 0;            | ne = 1; (state size)
%        x(s1=1) = 0;       x(s2=1) = 0;            |
%   Use states x1 = x,    x2 = x_{t}                |
%       =>                                          |
%   PDE: x1_{t} = x2;                               | 
%        x2_{t} = c1*x1_{(2,0)} + c2*x1_{(0,2)};    |
%   BCs: x1(s1=0) = 0;      x1(s2=0) = 0;           |
%        x1(s1=1) = 0;      x1(s2=1) = 0;           |
%        x2(s1=0) = 0;      x2(s2=0) = 0;           |
%        x2(s1=1) = 0;      x2(s2=1) = 0;           |
pvar s1 s2 t; % define independent variables
pdevar x1 x2; x1 = setInfinite(x1,2); x2 = setInfinite(x2,2);
c1 = 1;   c2 = 1;
%% Define equations
%   PDE: x1_{t} = x2;                               | 
%        x2_{t} = c1*x1_{(2,0)} + c2*x1_{(0,2)};    |
eq_PDE = [diff(x1,t)==x2;
                 diff(x2,t)==c1*diff(x1,s1,2)+c2*diff(x1,s2,2)]; 
%   BCs: x1(s1=0) = 0;      x1(s2=0) = 0;           |
%        x1(s1=1) = 0;      x1(s2=1) = 0;           |
%        x2(s1=0) = 0;      x2(s2=0) = 0;           |
%        x2(s1=1) = 0;      x2(s2=1) = 0;           | 
eq_BC = [subs(x1,s1,0)==0;subs(x1,s1,1)==0;subs(x1,s2,0)==0;subs(x1,s2,1)==0;
            subs(x2,s1,0)==0;subs(x2,s1,1)==0;subs(x2,s2,0)==0;subs(x2,s2,1)==0]; 
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde',[eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
