clc; clear; clear stateNameGenerator
pvar s1 s2 t; % define independent variables

%% Define dependent variables and system variable
%   Telegraph equation:
%   PDE: x1_{t} = x2
%        x2_{t} = -2*lam*x2 + C^2*(x1_{s1s1} + x1_{s2s2})
%   BCs: x1 = 0 and x2 = 0 on s1=0, s1=1, s2=0, and s2=1
%   C = 1; lam = 1
C = 1;
lam = 1;
x1 = pde_var(1,[s1;s2],[0,1;0,1]);
x2 = pde_var(1,[s1;s2],[0,1;0,1],[2;2]);

%% Define equations
eq_PDE = [diff(x1,'t')==x2;
          diff(x2,'t')==-2*lam*x2+C^2*(diff(x1,s1,2)+diff(x1,s2,2))];
eq_BC = [subs([x1;x2],s1,0)==0;
         subs([x1;x2],s1,1)==0;
         subs([x1;x2],s2,0)==0;
         subs([x1;x2],s2,1)==0];

%% initialize pde system;
PDE = initialize([eq_PDE; eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);
