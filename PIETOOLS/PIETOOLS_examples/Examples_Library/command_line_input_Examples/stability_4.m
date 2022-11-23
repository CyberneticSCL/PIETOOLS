clc; clear;
pvar s t theta; % define independent variables
%% Define dependent variables and system variable
% 	PDE: x1_{t} = sig1*x2 - (1/r1)*x1_{s}           | Different parameters          Saba 2019 [3] 
%        x2_{t} = sig2*x1 + (1/r2)*x2_{s}           | may be invoked calling  
%   BCs: x1(s=0) = qb*x2(s=0)                       | examples 4.1, 4.2, 
%        x2(s=1) = pb*x1(s=1)                       | 4.3, 4.4, and 4.5.
x = state('pde',2); % x = [x1; x2];
pde = sys();
r1=.8;   r2=1.1;   sig1=2.3;   sig2=-3.5;    qb=-.7;   pb=.5; 
%% Define equations
eq_PDE = diff(x,t)==[0,sig1;sig2,0]*x+[-1/r1,0;0,1/r2]*diff(x,s); % 	PDE: x1_{t} = sig1*x2 - (1/r1)*x1_{s},x2_{t} = sig2*x1 + (1/r2)*x2_{s}
eq_BC = [[1,0]*subs(x,s,0); [0,1]*subs(x,s,1)]==[0,qb;pb,0]*[[1,0]*subs(x,s,0);[0,1]*subs(x,s,1)];    
%   BCs: x1(s=0) = qb*x2(s=0), x2(s=1) = pb*x1(s=1)                      
%% addequations to pde system; set control inputs/observed inputs, if any
pde = addequation(pde,[eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
