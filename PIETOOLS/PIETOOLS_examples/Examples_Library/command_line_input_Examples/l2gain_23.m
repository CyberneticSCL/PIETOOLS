%clc; clear; clear stateNameGenerator
pvar s t; % define independent variables

%% Define dependent variables and system variable
%   PDE: x_{t} = A2(s)*x_{ss} + A1(s)*x_{s} + A0(s)*x + w(t)
%   BCs: x(s=0) = 0,        x_{s}(s=1) = 0
%   Out: z(t) = x(t,1)
%   A2 = s^3 - s^2 + 2;   A1 = 3*s^2 - 2*s
%   A0 = -0.5*s^3 + 1.3*s^2 - 1.5*s + 0.7 + lam;  lam = 4.6
pde_var x output z input w
x.vars = s;     x.dom = [0,1];
lam = 4.6;
A2 = s^3-s^2+2;
A1 = 3*s^2-2*s;
A0 = -0.5*s^3+1.3*s^2-1.5*s+0.7+lam;

%% Define equations
eq_PDE = diff(x,t)==A0*x+A1*diff(x,s)+A2*diff(x,s,2)+w;
eq_BC = [subs(x,s,0)==0; subs(diff(x,s),s,1)==0];
eq_out = z==subs(x,s,1);

%% initialize pde system;
PDE = initialize([eq_PDE; eq_BC; eq_out]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);
