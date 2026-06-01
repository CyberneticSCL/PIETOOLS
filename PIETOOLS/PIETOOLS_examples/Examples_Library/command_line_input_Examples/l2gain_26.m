clc; clear; clear stateNameGenerator
pvar s t; % define independent variables

%% Define dependent variables and system variable
%   PDE: x_{t} = Cm(s)*x + (1/R)*x_{ss} + [s;s;s]*w(t)
%   BCs: x(s=0) = 0,        x(s=1) = 0
%   Out: z(t) = int(x(t,s),s,0,1)
%   Cm = [0,0,0; s,0,0; s^2,-s^3,0];   R = 21 - 1e-3
pde_var x output z input w
x.vars = s;     x.dom = [0,1];     x.size = 3;
z.size = 3;
Cm = [0,0,0; s,0,0; s^2,-s^3,0];
R = 21-1e-3;

%% Define equations
eq_PDE = diff(x,t)==Cm*x+(1/R)*diff(x,s,2)+s*ones(3,1)*w;
eq_BC = [subs(x,s,0)==0; subs(x,s,1)==0];
eq_out = z==int(x,s,[0,1]);

%% initialize pde system;
PDE = initialize([eq_PDE; eq_BC; eq_out]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);
