%   ODE: xo_{t} = xo + x_{s}(s=0)               	| k = -2                (stable for k=-2)            
%   PDE: x_{t} = x_{ss}                             |                               Tang 2011 [13]
%   BCs: x(s=0) = xo,    x(s=1) = k*xo              |                           (PDE stable )
clc; clear; clear stateNameGenerator
pvar s t; % define independent variables
% Define dependent variables and system variable
x=pde_var(); X = pde_var(s,[0,1]);
%% Specify the parameters
k = -2;
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%% Define equations
eq_ODE=diff(x,t)==x+subs(diff(X,s),s,0); % ODE: xo_{t} = xo + x_{s}(s=0)
eq_PDE = diff(X,t)==diff(X,s,2); % 	PDE: x_{t} = x_{ss}
eq_BC = [subs(X,s,0)==-x;subs(X,s,1)==k*x]; %  BCs: x(s=0) = -xo,    x(s=1) = k*xo
%% initialize pde system;
PDE = initialize([eq_ODE;eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);



% %% Define dependent variables and system variable
% %   ODE: xo_{t} = xo + x_{s}(s=0)               	| k = -2   
% %   PDE: x_{t} = x_{ss}                             |                            
% %   BCs: x(s=0) = -xo,    x(s=1) = k*xo             |
% x=state('ode'); X = state('pde');
% pde = sys(); 
% k = -2;
% %% Define equations
% eq_ODE=diff(x,t)==x+subs(diff(X,s),s,0); % ODE: xo_{t} = xo + x_{s}(s=0)
% eq_PDE = diff(X,t)==diff(X,s,2); % 	PDE: x_{t} = x_{ss}
% eq_BC = [subs(X,s,0)==-x;subs(X,s,1)==k*x]; %  BCs: x(s=0) = -xo,    x(s=1) = k*xo
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_ODE;eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
