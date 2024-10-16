clc; clear;
%% Define dependent variables and system variable
%   PDE: x_{t} = Cm(s)*x + (1/R)*x_{ss} + s*w(t)    | R = (21-1e-3)         (gamma =  4.23, for R = 21-1e-3)
%   BCs: x(s=0) = 0,        x(s=1) = 0              | Cm = [0,0,0;                  Shivakumar 2019 [12] (Example 2)
%   Out: z(t) = int(x(t,s),s,0,1)                   |       s,0,0;
%                                                   |       s^2,-s^3,0]
pvar s t; % define independent variables
outputvar z; inputvar w; pdevar x;
x.len = 3; z.len = 3;
R = 20;Cm = [0,0,0; s,0,0; s^2, -s^3,0];
%% Define equations
%   PDE: x_{t} = Cm(s)*x + (1/R)*x_{ss} + s*w(t) 
eq_PDE = diff(x,t)==Cm*x+(1/R)*diff(x,s,2)+s*[1;1;1]*w;
eq_BC = [subs(x,s,0)==0;subs(x,s,1)==0]; %  BCs: x(s=0) = 0,     x(s=1) = 0 
eq_out = [int(x,s,[0,1])];    %   Out: z(t) = int(x(t,s),s,0,1) 
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde',[eq_PDE;eq_BC;eq_out]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
