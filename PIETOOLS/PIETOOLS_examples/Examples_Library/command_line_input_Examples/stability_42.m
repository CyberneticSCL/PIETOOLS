clc; clear; clear stateNameGenerator
pvar s1 s2 t; % define independent variables

%% Define dependent variables and system variable
%   ODE: x1_{t} = A*x1 + A1*x3(t,1) + B*x2(t,s1=0)
%   PDE: x2_{t} = x2_{s1s1} + a*x2 + a2*x4(t,s1,s2=1)
%        x3_{t} = x3_{s2}
%        x4_{t} = x4_{s2}
%   tau = 1; A = -1; A1 = 0.5; B = -1; a = 1; a2 = 0.5
tau = 1;
A = -1;
A1 = 0.5;
B = -1;
a = 1;
a2 = 0.5;
pde_var x1 x2 x3 x4
x2.vars = s1;       x2.dom = [0,1];
x3.vars = s2;       x3.dom = [1,1+tau];
x4.vars = [s1;s2];  x4.dom = [0,1;1,1+tau];

%% Define equations
eq_PDE = [diff(x1,'t')==A*x1+B*subs(x2,s1,0)+A1*subs(x3,s2,1);
          diff(x2,'t')==diff(x2,s1,2)+a*x2+a2*subs(x4,s2,1);
          diff(x3,'t')==diff(x3,s2);
          diff(x4,'t')==diff(x4,s2)];
eq_BC = [subs(diff(x2,s1),s1,0)==0;
         subs(x2,s1,1)==0;
         subs(x3,s2,1+tau)==x1;
         subs(diff(x4,s1),s1,0)==0;
         subs(x4,s1,1)==0;
         subs(x4,s2,1+tau)==x2];

%% initialize pde system;
PDE = initialize([eq_PDE; eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);
