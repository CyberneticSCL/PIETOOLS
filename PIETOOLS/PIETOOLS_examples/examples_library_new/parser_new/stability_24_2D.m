clc; clear;
%% Define dependent variables and system variable
%   PDE: x_{t} = a*x                                | a = 4.0;              (stable for a <= 0.5*pi^2) 
%                 + b1*x_{(1,0)} + b2*x_{(0,1)}     | b1 = 0;   b2 = 0;             Demetriou, 2019 [17]   
%                  + c1*x_{(2,0)} + c2*x_{(0,2)}    | c1 = 1;   c2 = 1;
%   BCs: x(s1=0) = 0,       x_{s1}(s1=1) = 0,       | ne = 1 (state size)
%        x(s2=0) = 0,       x_{s2}(s2=1) = 0;       |
pvar s1 s2 t; % define independent variables
pdevar x1; x1 = setInfinite(x1,2); 
a=4;c1 = 1; c2 =1; b1 = 0; b2 = 0;
%% Define equations
%   PDE: x_{t} = a*x                                | a = 4.0;              (stable for a <= 0.5*pi^2) 
%                 + b1*x_{(1,0)} + b2*x_{(0,1)}     | b1 = 0;   b2 = 0;             Demetriou, 2019 [17]   
%                  + c1*x_{(2,0)} + c2*x_{(0,2)}    | c1 = 1;   c2 = 1;
eq_PDE = [diff(x1,t)==a*x1+c1*diff(x1,s1,2)+c2*diff(x2,s2,2)+b1*diff(x1,s1)+b2*diff(x1,s2)]; 
%   BCs: x1(s1=0) = 0,   x1(s2=0) = 0,              | 
%        x1(s1=1) = 0,   x1(s2=1) = 0;              |

eq_BC = [subs(x1,s1,0)==0;subs(x1,s1,1)==0;
                 subs(x1,s2,0)==0;subs(x1,s2,1)==0;]; 
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde',[eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
