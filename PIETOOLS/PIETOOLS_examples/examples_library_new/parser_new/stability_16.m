clc; clear;
%% Define dependent variables and system variable
%   PDE: 0 = alp*w_{ssss} + bet*w_{tt}              |
%             - gam*w_{ttss} + w_{tttt}             |
%   Use states:                                     |
%        x1 = w_{ttt},      x2 = w_{t},             |
%        x3 = w_{tt},       x4 = w.                 |
%       =>                                          |
%   PDE: x1_{t} = -alp*x4_{ssss} - bet*x3           |
%                 + gam*x3_{ss}                     |
%        x2_{t} = x3,  x3_{t} = x1,  x4_{t} = x2    |
%   BCs: w(s=0) = 0,        w_{s}(s=0) = 0          |
%        w_{ss}(s=1) - w(s=1) = 0                   |
%        w_{sss}(s=1) - w_{s}(s=1) = 0              |
%        w_{t}(s=0) = 0     w_{ts}(s=0) = 0         |
%        w_{tt}(s=0) = 0    w_{tts}(s=0) = 0        |
pvar s t; % define independent variables
pdevar x1 x2 x3 x4;
%% Define equations
%   PDE: x1_{t} = -alp*x4_{ssss} - bet*x3           |
%                 + gam*x3_{ss}                     |
%        x2_{t} = x3,  x3_{t} = x1,  x4_{t} = x2    |
eq_PDE = [diff(x1,t)==-alp*diff(x4,s,4)-bet*x3+gam*diff(x2,s,2); 
          diff(x2,t)==x3;
          diff(x3,t)==x1;
          diff(x4,t)==x2];
%   BCs: x4(s=0) = 0,        x4_{s}(s=0) = 0          |
%        x4_{ss}(s=1) - x4(s=1) = 0                   |
%        x4_{sss}(s=1) - x4_{s}(s=1) = 0              |
%        x3(s=0) = 0     x3_{s}(s=0) = 0         |
%        x2(s=0) = 0    x2_{s}(s=0) = 0        |
eq_BC = [subs(x4,s,0)==0;subs(x4,s,0)==0;
             subs(diff(x4,s,2),s,1)-subs(x4,s,1)==0;
             subs(diff(x4,s,3),s,1)-subs(diff(x4,s),s,1)==0;
             subs(x3,s,0)==0; subs(x3,s,0)==0; 
             subs(x2,s,0)==0; subs(diff(x2,s),s,0)==0]; 
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde',[eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
