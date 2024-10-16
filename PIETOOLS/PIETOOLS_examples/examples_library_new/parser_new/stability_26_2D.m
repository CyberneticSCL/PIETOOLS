clc; clear;
%% Define dependent variables and system variable
%   PDE: p_{t}  = -s2*p_{s1} - v1_{s1} - v2_{s2}    | M = 0.1;
%        v1_{t} = -s2*v1_{s1} -v2 -(1/M^2)*p_{s1}   | lam = 1; nu = 1;              Antonelli, 2021 [17]
%                  +nu*(v1_{s1s1} + v1_{s2s2})      |
%                   +lam*(v1_{s1s1} + v2_{s2s1})    |
%        v2_{t} = -s2*v2_{s1} -(1/M^2)*p_{s2}       |
%                  +nu*(v2_{s1s1} + v2_{s2s2})      |
%                   +lam*(v1_{s1s2} + v2_{s2s2})    |
%   BCs: p(s1=0) = 0;      p(s2=0) = 0;             |
%        v(s1=0) = 0;      v(s2=0) = 0;             |
%        v(s1=1) = 0;      v(s2=1) = 0;             |
pvar s1 s2 t; % define independent variables
pdevar p v1 v2; v1 = setInfinite(v1,2); v2 = setInfinite(v2,2); p = setInfinite(p,2);
lam = 1; nu = 1; M = 0.1;
%% Define equations
%   PDE: p_{t}  = -s2*p_{s1} - v1_{s1} - v2_{s2}    | M = 0.1;
%        v1_{t} = -s2*v1_{s1} -v2 -(1/M^2)*p_{s1}   | lam = 1; nu = 1;              Antonelli, 2021 [17]
%                  +nu*(v1_{s1s1} + v1_{s2s2})      |
%                   +lam*(v1_{s1s2} + v2_{s2s1})    |
%        v2_{t} = -s2*v2_{s1} -(1/M^2)*p_{s2}       |
%                  +nu*(v2_{s1s1} + v2_{s2s2})      |
%                   +lam*(v1_{s1s2} + v2_{s2s2})    |
eq_PDE = [diff(p,t)==-s2*diff(p,s1)-diff(v1,s1)-diff(v2,s2);
                 diff(v1,t)==-s2*diff(v1,s1)-v2-(1/M^2)*diff(p,s1)+nu*(diff(v1,s1,2) + diff(v1,s2,2))+lam*(diff(diff(v1,s1),s1) + diff(diff(v2,s2),s1));
                 diff(v2,t)==-s2*diff(v2,s1)-(1/M^2)*diff(p,s2)+nu*(diff(v2,s1,2) + diff(v2,s2,2))+lam*(diff(diff(v2,s1),s2) + diff(diff(v2,s2),s2))]; 
%   BCs: p(s1=0) = 0;      p(s2=0) = 0;             |
%        v(s1=0) = 0;      v(s2=0) = 0;             |
%        v(s1=1) = 0;      v(s2=1) = 0;             |
eq_BC = [subs([v1;v2],s1,0)==0;subs([v1;v2],s1,1)==0;subs([v1;v2],s2,0)==0;subs([v1;v2],s2,1)==0;
            subs(p,s1,0)==0;subs(p,s2,0)==0]; 
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde',[eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
