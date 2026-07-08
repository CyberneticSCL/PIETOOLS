% Tip-damped wave eq with reaction term and viscous damping
%     | PDE: u_{tt} + 2*ad*u_{t} = -ad^2*u + u_{ss} | k = 1     Datko 1986 [9] (Test 7.5c)
%     | BCs: u(s=0) = 0                             |
%     |      u_{s}(s=1) = -k*u_{t}(s=1)             |
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
%   Use states x1 = u_{t}, x2 = u.                  | k = 1         (PDE exponentially dual stable )                 
%       =>                                          | ad = 1                    (Finite Energy PIE to PDE stable )   
%   PDE: x1_{t} = -2*ad*x1 - ad^2*x2 + x2_{ss}      |
%   BCs: x1(0) = 0,     x2(0) = 0,                  |
%        k*x1(1) + x2_{s}(1) = 0                    |
%clc; clear;
clear stateNameGenerator
pvar s t; % define independent variables
%% Specify parameters
k = 1;    ad = 1;
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%% Define dependent variables and system variable
pde_var x1 x2
x1.vars = s;    x2.vars = s; 
%% Define equations
eq_PDE = [diff(x1,t,1)==-2*ad*x1-ad^2*x2+diff(x2,s,2); %   PDE: x1_{t} = -2*ad*x1 - ad^2*x2 + x2_{ss}      |
                  diff(x2,t,1)==x1;];% x2_{t} = x1 
eq_BC = [subs(x1,s,0)==0;subs(x2,s,0)==0;%x1(0) = 0,     x2(0) = 0,                  |
                k*subs(x1,s,1)+subs(diff(x2,s,1),s,1)==0;]; %        k*x1(1) + x2_{s}(1) = 0      
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);