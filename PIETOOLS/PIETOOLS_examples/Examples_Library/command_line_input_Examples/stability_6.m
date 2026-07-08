%   | PDE: x_{t} = lam*x + x_{ss}                   | lam = 2.466           (unstable for lam > pi^2/4 \approx 2.467,  
%                                                                                                     exponentially stable and FE PDE stable for lam < pi^2/4)
%   | BCs: x(s=0) = 0,      x_{s}(s=1) = 0          |                               Gahlawat 2017 [4]
clc; clear; clear stateNameGenerator
pvar s t; % define independent variables
% Define dependent variables and system variable
x = pde_var(s,[0,1]);
%% Specify the parameters
lam = 2.466;
 if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
 end
%% Define equations
eq_PDE = diff(x,t)==lam*x+diff(x,s,2); % 	PDE: x_{t} = lam*x + x_{ss}
eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,1)==0]; %  BCs: x(s=0) = 0,  x_{s}(s=1) = 0
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE2 = convert(PDE,'pie');
display(PIE2);



% %% Define dependent variables and system variable
% %   | PDE: x_{t} = lam*x + x_{ss}                   | lam = 2.466 
% %   | BCs: x(s=0) = 0,      x_{s}(s=1) = 0
% x = state('pde');
% pde = sys(); lam = 2.466;
% %% Define equations
% eq_PDE = diff(x,t)==lam*x+diff(x,s,2); % 	PDE: x_{t} = lam*x + x_{ss}
% eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,1)==0]; %  BCs: x(s=0) = 0,  x_{s}(s=1) = 0
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
