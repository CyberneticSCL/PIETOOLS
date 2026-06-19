%% Euler-Bernoulli Equation
% using PDE states x1=w_{t}, x2=w and PIE states (w_{t},w_{ssss})
%clc; clear; 
clear stateNameGenerator
pvar s t; % define independent variables
% TImoshenko Beam Model, 4-th order in time implementation
%   PDE: 0 = alp*w_{ssss} + bet*w_{tt}              |
%             - gam*w_{ttss} + delt*w_{tttt}             |
% alp=E*I, beta= rho*A, gam= rhoI+E*I*rho/k/G, rho^2*I/k/G 
%   Use states:                                     |
%        x1 = w_{ttt},      x2 = w_{t},             |
%        x3 = w_{tt},       x4 = w.                 |
% Specify the parameters
% Assuming E=I=rho=A=G=k=1,
alp = 1;     bet = 1;
% Case gam=delt=0 collapse to Euler-Bernoulli model,
% with PDE states x1=w_{t}, x2=w and PIE states (w_{t},w_{ssss})
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
%   We rewrite the system as a single equation:     |                     
%   PDE: 0 = alp*w_{ssss} + bet*w_{tt}-c*w(t)      |
%     BCs w(s=0) = 0,        w_{ss}(s=1) = 0         |                               Peet 2019 [8] (Example 8.1.0.1)
%        w_{s}(s=0) = 0,    w_{sss}(s=1) = 0        |
%   Use states:                                     |
%        x1 = w_{t},      x2 = w.                 |
%       =>                                       
%   PDE: x1_{t} = -alp/bet*x2_{ssss},
%             x2_{t} = x_1           |
%   BCs: x2(s=0) = 0,        x2_{ss}(s=1) = 0      |
%        x2_{s}(s=0)= 0        x2_{sss}(s=1)= 0    |
pde_var x1 x2 %x2 x3 x4
c=0.1;
x1.vars = s;   x2.vars = s;   % x3.vars = s;    x4.vars = s;   
%% Define equations
eq_PDE = [diff(x1,t,1)==-alp/bet*diff(x2,s,4)-c*x1; %  %   PDE: x1_{t} = -alp/bet*x2_{ssss}-c/bet*x1,
                    diff(x2,t,1)==x1;]; %             x2_{t} = x_1           |
eq_BC = [subs(x2,s,0)==0;
   % subs(x1,s,0)==0; 
    subs(diff(x2,s,2),s,1)==0;%%   BCs: x2(s=0) = 0,        x2_{ss}(s=1) = 0      |
                subs(diff(x2,s,1),s,0)==0; 
                subs(diff(x2,s,3),s,1)==0;];%        x2_{s}(s=0)= 0        x2_{sss}(s=1)= 0    |
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);
