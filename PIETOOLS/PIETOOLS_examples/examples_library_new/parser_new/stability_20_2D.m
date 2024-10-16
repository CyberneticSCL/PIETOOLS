clc; clear;
%% Define dependent variables and system variable
%   PDE: x_{t} = c1*x_{s1} + c2*x_{s2}              | c1 = 1; c2 = 1;
%   BCs: x(s1=0) = 0,   x(s2=0) = 0;                |   ne = 1 (state size)
pvar s1 s2 t; % define independent variables
pdevar x; x = setInfinite(x,2);
c1 = 1; c2 = 1;
%% Define equations
eq_PDE = diff(x,t)==c1*diff(x,s1)+c2*diff(x,s2) ; %   PDE: x_{t} = c1*x_{s1} + c2*x_{s2}  
eq_BC = [subs(x,s1,0)==0;subs(x,s2,0)==0];        %   BCs: x(s1=0) = 0,   x(s2=0) = 0;                
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde', [eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
