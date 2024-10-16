clc; clear;
%% Define dependent variables and system variable
%   PDE: x_{t} = C*(x_{(2,0} + x_{(0,2)})           | C = 1;
%                 - b1*x_{(1,0)} - b2*x_{(0,1)}     | b1 = 0.5; b2 = 2;             Holmes, 1994 [14] Eq. (2)
%   BCs: x(s1=0) = 0,   x(s2=0) = 0,                | ne = 1 (state size)
%        x(s1=1) = 0,   x(s2=1) = 0;                |
pvar s1 s2 t; % define independent variables
pdevar x; x = setInfinite(x,2);
b1 = 0.5; b2 = 2;C = 1;
%% Define equations
%   PDE: x_{t} = C*(x_{(2,0} + x_{(0,2)})           | C = 1;
%                 - b1*x_{(1,0)} - b2*x_{(0,1)}     | b1 = 0.5; b2 = 2;             Holmes, 1994 [14] Eq. (2)
eq_PDE = diff(x,t)==C*diff(x,s1,2)+C*diff(x,s2,2)-b1*diff(x,s1)-b2*diff(x,s2); 
%   BCs: x(s1=0) = 0,   x(s2=0) = 0,                | ne = 1 (state size)
%        x(s1=1) = 0,   x(s2=1) = 0;  
eq_BC = [subs(x,s1,0)==0;subs(x,s1,1)==0;subs(x,s2,0)==0;subs(x,s2,1)==0]; 
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde',[eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
