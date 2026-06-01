clc; clear; clear stateNameGenerator
pvar s1 s2 t; % define independent variables

%% Define dependent variables and system variable
%   Linearized 2D isentropic compressible Navier-Stokes:
%   p_{t}  = -s2*p_{s1} - v1_{s1} - v2_{s2}
%   v1_{t} = -s2*v1_{s1} - v2 - (1/M^2)*p_{s1}
%            + nu*(v1_{s1s1}+v1_{s2s2}) + lam*(v1_{s1s1}+v2_{s2s1})
%   v2_{t} = -s2*v2_{s1} - (1/M^2)*p_{s2}
%            + nu*(v2_{s1s1}+v2_{s2s2}) + lam*(v1_{s1s2}+v2_{s2s2})
%   BCs: p(s1=0)=0, p(s2=0)=0, and v=0 on all edges
M = 0.1;
lam = 1;
nu = 1;
v1 = pde_var([s1;s2],[0,1;0,1]);
v2 = pde_var([s1;s2],[0,1;0,1]);
p = pde_var([s1;s2],[0,1;0,1]);

%% Define equations
eq_PDE = [diff(p,'t')==-s2*diff(p,s1)-diff(v1,s1)-diff(v2,s2);
          diff(v1,'t')==-s2*diff(v1,s1)-v2-(1/M^2)*diff(p,s1) ...
              +nu*(diff(v1,s1,2)+diff(v1,s2,2)) ...
              +lam*(diff(v1,s1,2)+diff(v2,[s1;s2]));
          diff(v2,'t')==-s2*diff(v2,s1)-(1/M^2)*diff(p,s2) ...
              +nu*(diff(v2,s1,2)+diff(v2,s2,2)) ...
              +lam*(diff(v1,[s1;s2])+diff(v2,s2,2))];
eq_BC = [subs([v1;v2],s1,0)==0;
         subs([v1;v2],s1,1)==0;
         subs([v1;v2],s2,0)==0;
         subs([v1;v2],s2,1)==0;
         subs(p,s1,0)==0;
         subs(p,s2,0)==0];

%% initialize pde system;
PDE = initialize([eq_PDE; eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);
