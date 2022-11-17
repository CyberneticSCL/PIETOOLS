clc; clear;
pvar s t theta; % define independent variables
%% Define dependent variables and system variable
%	PDE: x_{t} = Fm*x - Lm*x_{s}                    | Different parameters          Lamare 2016 [1] 
%   BCs: [x_-(s=1)] = [Gm1, Gm2] [x_-(s=0)]         | may be invoked calling
%        [x_+(s=0)] = [Gm3, Gm4] [x_+(s=1)]         | examples 2.1, 2.2, 
%                                                   | 2.3, and 2.4.
x = state('pde',2); % x = [x_-; x_+];
pde = sys();
Gm1=[.2]; Gm2=[-.3]; Gm3=[.6]; Gm4=[.1]; Lm=[-3 0;0 1]; Fm=[.2 -.3; .6 .1];
%% Define equations
eq_PDE = diff(x,t)==Fm*x-Lm*diff(x,s); % 	PDE: x_{t} = Fm*x - Lm*x_{s}
eq_BC = [[1,0]*subs(x,s,1); [0,1]*subs(x,s,0)]==[Gm1,Gm2;Gm3,Gm4]*[[1,0]*subs(x,s,0);[0,1]*subs(x,s,1)];    
%   BCs: [x_-(s=1)] = [Gm1, Gm2] [x_-(s=0)], [x_+(s=0)] = [Gm3, Gm4] [x_+(s=1)]
%% addequations to pde system; set control inputs/observed inputs, if any
pde = addequation(pde,[eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
