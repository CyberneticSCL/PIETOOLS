clc; clear; clear stateNameGenerator
pvar s t; % define independent variables
%
%   PDE: r*aa * w_{tt} = k*aa*g * (-phi_{s} + w_{ss})                               Peet 2019 [8] (Example 8.1.0.2)
%        r*II * phi_{tt} = E*II * phi_{ss}  + k*aa*g * (w_{s} - phi) 
%   BCs: phi(s=0) = 0,      phi_{s}(s=1) = 0    
%        w(s=0) = 0,        w_{s}(s=1) - phi(s=1) = 0           
% 
%% Specify the parameters
alp = 1;     bet = 1;     gam = 2;
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
%   We rewrite the system as a single equation:     |                       4th order implementation (stable)
%   PDE: 0 = alp*w_{ssss} + bet*w_{tt}              |
%             - gam*w_{ttss} + w_{tttt}             |
%   Use states:                                     |
%        x1 = w_{ttt},      x2 = w_{t},             |
%        x3 = w_{tt},       x4 = w.                 |
%       =>                                        
%   PDE: x1_{t} = -alp*x4_{ssss} - bet*x3           |
%                 + gam*x3_{ss}                     |
%        x2_{t} = x3,  x3_{t} = x1,  x4_{t} = x2    |
%   BCs: x4(s=0) = 0,        x4_{s}(s=0) = 0          |
%        x4_{ss}(s=1) - x4(s=1) = 0                   |
%        x4_{sss}(s=1) - x4_{s}(s=1) = 0              |
%        x2(s=0) = 0     x2_{s}(s=0) = 0         |
%        x3(s=0) = 0    x3_{s}(s=0) = 0        |
%% Define dependent variables and system variable
pde_var x1 x2 x3 x4
x1.vars = s;    x2.vars = s;    x3.vars = s;    x4.vars = s;   
%% Define equations
eq_PDE = [diff(x1,t,1)==-alp*diff(x4,s,4)-bet*x3+gam*diff(x3,s,2); %   PDE: x1_{t} = -alp*x4_{ssss} - bet*x3   + gam*x3_{ss}         |
                  diff(x2,t,1)==x3; %    x2_{t} = x3,  
                  diff(x3,t,1)==x1;% x3_{t} = x1,  
                  diff(x4,t,1)==x2;];% x4_{t} = x2    |
eq_BC = [ subs(x1,s,0)==0; subs(x2,s,0)==0;
                subs(x3,s,0)==0;subs(x4,s,0)==0; 
                subs(diff(x4,s,1),s,0)+int(x3,s,0,1)==0;%  x4(s=0) = 0,        x4_{s}(s=0) = 0          |
                subs(diff(x4,s,3),s,0)-subs(diff(x3,s,1),s,0)+subs(diff(x4,s,1),s,0)==0; % x4_{ss}(s=1) - x4(s=1) = 0                   |
                subs(diff(x4,s,2),s,1)- subs(x3,s,1)==0;   %        x4_{sss}(s=1) - x4_{s}(s=1) = 0              |
                subs(diff(x4,s,3),s,1)-gam*subs(diff(x3,s,1),s,1)] %        x3(s=0) = 0    x3_{s}(s=0) = 0        |
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);
