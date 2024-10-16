clc; clear;
%% Define dependent variables and system variable
%   ODE: x1_{t} = (A+BK)*x1 + B*x2(s1=0,s2=0)       | A = I; B = I;
%   PDE: x2_{t} = c1*x_{(2,0)} + c2*x_{(0,2)}       | K = -2*I;
%   BCs: x2_{(1,0)}(s1=0) = 0;  x2(s1=1) = 0;       |
%        x2_{(0,1)}(s2=0) = 0;  x2(s2=1) = 0;       |
pvar s1 s2 t; % define independent variables
odevar x1; pdevar x2; x2 = setInfinite(x2,2); 
I = eye(1); A = I; B = I; K = -2*I;
%% Define equations
%   ODE: x1_{t} = (A+BK)*x1 + B*x2(s1=0,s2=0)       | A = I; B = I;
%   PDE: x2_{t} = c1*x_{(2,0)} + c2*x_{(0,2)}       | K = -2*I;
eq_ODE = diff(x1,t)==(A+BK)*x1+B*subs(subs(x2,s1=0),s2=0);
eq_PDE = [diff(x2,t)==c1*diff(x2,s1,2)+c2*diff(x2,s2,2)]; 
%   BCs: x2_{(1,0)}(s1=0) = 0;  x2(s1=1) = 0;       |
%        x2_{(0,1)}(s2=0) = 0;  x2(s2=1) = 0;       |
eq_BC = [subs(diff(x2,s1),s1,0)==0;subs(x2,s1,1)==0;
                subs(diff(x2,s2),s2,0)==0;subs(x2,s2,1)==0]; 
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde',[eq_ODE;eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
