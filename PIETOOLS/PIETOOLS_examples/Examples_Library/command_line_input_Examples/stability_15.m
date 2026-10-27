% Timoshenko Beam
%   PDE: r*aa * w_{tt} = k*aa*g * (-phi_{s} + w_{ss})                               Peet 2019 [8] (Example 8.1.0.2)
%        r*II * phi_{tt} = E*II * phi_{ss}  + k*aa*g * (w_{s} - phi) 
%   BCs: phi(s=0) = 0,      phi_{s}(s=1) = 0    
%        w(s=0) = 0,        w_{s}(s=1) - phi(s=1) = 0           
%   Assume all parameters are 1, and use states:    | Hyperbolic/diffusive implementation (unstable)
%        x1 = w_{t},    x2 = w_{s},                 |
%        x3 = phi_{t},  x4 = phi.                   |                  (PDE stable with damping c >0)
%       =>                                          |
%   PDE: x1_{t} = x2_{s} - x4_{s}                   |
%        x2_{t} = x1_{s}                            |
%        x3_{t} = x2 - x4 -c x3                         |
%        x4_{t} = x3                                |
%   BCs: x4(0) = 0,     x4_{s}(1) = 0,              |
%        x3(0) = 0,     x1(0) = 0,                  |
%        x2(1) - x4(1) = 0 
clc; clear; clear stateNameGenerator
pvar s t; % define independent variables
%% Specify the parameters 
k = 1;    aa = 1;    II = 1;    g = 1;   
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%% Define dependent variables and system variable
% use PDE states:  (unstable Hyperbolic/diffusive implementation)
%        x1 = w_{t},    x2 = w_{s},                 |
%        x3 = phi_{t},  x4 = phi.                   |
%   PDE: x1_{t} = x2_{s} - x4_{s}                   |
%        x2_{t} = x1_{s}                            |
%        x3_{t} = x2 - x4                           |
%        x4_{t} = x3                                |
%   BCs: x4(0) = 0,     x4_{s}(1) = 0,              |
%        x3(0) = 0,     x1(0) = 0,                  |
%        x2(1) - x4(1) = 0           
pde_var x1 x2 x3 x4
x1.vars = s;    x2.vars = s;    x3.vars = s;    x4.vars = s;   
%% Define equations
%   PDE: x1_{t} = x2_{s} - x4_{s}                   |
%        x2_{t} = x1_{s}                            |
%        x3_{t} = x2 - x4                           |
%        x4_{t} = x3                                |
eq_PDE = [diff(x1,t)==diff(x2,s)-diff(x4,s); 
          diff(x2,t)==diff(x1,s);
          diff(x3,t)==x2-x4;
          diff(x4,t)==x3];
%   BCs: x4(0) = 0,     x4_{s}(1) = 0,              |
%        x3(0) = 0,     x1(0) = 0,                  |
%        x2(1) - x4(1) = 0 
eq_BC = [subs(x4,s,0)==0;subs(x3,s,0)==0;subs(diff(x4,s),s,1)==0;
         subs(x1,s,0)==0;subs(x2,s,1)-subs(x4,s,1)==0]; 
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);



% %% Define dependent variables and system variable
% %   PDE: x1_{t} = x2_{s} - x4_{s}                   |
% %        x2_{t} = x1_{s}                            |
% %        x3_{t} = x2 - x4                           |
% %        x4_{t} = x3                                |
% %   BCs: x4(0) = 0,     x4_{s}(1) = 0,              |
% %        x3(0) = 0,     x1(0) = 0,                  |
% %        x2(1) - x4(1) = 0           
% x1=state('pde');x2=state('pde');x3=state('pde');x4=state('pde'); 
% pde = sys();     
% %% Define equations
% %   PDE: x1_{t} = x2_{s} - x4_{s}                   |
% %        x2_{t} = x1_{s}                            |
% %        x3_{t} = x2 - x4                           |
% %        x4_{t} = x3                                |
% eq_PDE = [diff(x1,t)==diff(x2,s)-diff(x4,s); 
%           diff(x2,t)==diff(x1,s);
%           diff(x3,t)==x2-x4;
%           diff(x4,t)==x3];
% %   BCs: x4(0) = 0,     x4_{s}(1) = 0,              |
% %        x3(0) = 0,     x1(0) = 0,                  |
% %        x2(1) - x4(1) = 0 
% eq_BC = [subs(x4,s,0)==0;subs(x3,s,0)==0;subs(diff(x4,s),s,1)==0;
%          subs(x1,s,0)==0;subs(x2,s,1)-subs(x4,s,1)==0]; 
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
