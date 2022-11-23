clc; clear;
pvar s t theta; % define independent variables
%% Define dependent variables and system variable
%   PDE: x1_{t} = x2_{s} - x4_{s}                   |
%        x2_{t} = x1_{s}                            |
%        x3_{t} = x2 - x4                           |
%        x4_{t} = x3                                |
%   BCs: x4(0) = 0,     x4_{s}(1) = 0,              |
%        x3(0) = 0,     x1(0) = 0,                  |
%        x2(1) - x4(1) = 0           
x1=state('pde');x2=state('pde');x3=state('pde');x4=state('pde'); 
pde = sys();     
%% Define equations
%   PDE: x1_{t} = x2_{s} - x4_{s}                   |
%        x2_{t} = x1_{s}                            |
%        x3_{t} = x2 - x4                           |
%        x4_{t} = x3                                |
eq_PDE = [diff(x1,t)==diff(x2,s)-diff(x4,s); 
          diff(x2,t)==diff(x1,s);
          diff(x3,t)==x2-x4;
          diff(x4,t)==x3];
%   BCs: x4(0) = 0,     x4_{s}(1) = 0,              |
%        x3(0) = 0,     x1(0) = 0,                  |
%        x2(1) - x4(1) = 0 
eq_BC = [subs(x4,s,0)==0;subs(x3,s,0)==0;subs(diff(x4,s),s,1)==0;
         subs(x1,s,0)==0;subs(x2,s,1)-subs(x4,s,1)==0]; 
%% addequations to pde system; set control inputs/observed inputs, if any
pde = addequation(pde,[eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
