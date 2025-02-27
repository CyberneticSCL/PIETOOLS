clc; clear; clear stateNameGenerator
pvar s t; % define independent variables

%% Define dependent variables and system variable
%   PDE: x_{t} = -x_{s} + w(t)                      |
%   BCs: x(s=0) = 0                                 |
%   Out: z(t) = int(x(t,s),s,0,1)                   |          
pde_var x output z input w
x.vars = s;     x.dom = [0,1];    
%% Define equations
eq_PDE = [diff(x,t)==-diff(x,s)+w]; %   PDE: x_{t} = -x_{s} + w(t)
eq_BC = [subs(x,s,0)==0]; %   BCs: x(s=0) = 0 
eq_out = [z==int(x,s,[0,1])]; %   Out: z(t) = int(x(t,s),s,0,1)
%% initialize pde system;
pde = initialize([eq_PDE;eq_BC;eq_out]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);



% %% Define dependent variables and system variable
% %   PDE: x_{t} = -x_{s} + w(t)                      |
% %   BCs: x(s=0) = 0                                 |
% %   Out: z(t) = int(x(t,s),s,0,1)                   |          
% x=state('pde'); z = state('out'); w=state('in');
% pde = sys();     
% %% Define equations
% eq_PDE = [diff(x,t)==-diff(x,s)+w]; %   PDE: x_{t} = -x_{s} + w(t)
% eq_BC = [subs(x,s,0)==0]; %   BCs: x(s=0) = 0 
% eq_out = [z==int(x,s,[0,1])]; %   Out: z(t) = int(x(t,s),s,0,1)
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC;eq_out]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
