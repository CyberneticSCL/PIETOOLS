%   PDE: x_{t} = lam*x + x_{ss}                     | lam = 9.86   (unstable for lam > pi^2/4,
%                                                                                          exponentially stable and FE PDE stable for lam < pi^2 = 9.8696)  
%   BCs: x(s=0) = 0,      x(s=1) = 0                |                               Ahmadi 2015 [5] 
%                                                                         (Exponentially PDE stable )
%                                                                         (Finite Energy FEPIE2PDE stable ) 
clc; clear; clear stateNameGenerator
pvar s t; % define independent variables
% Define dependent variables and system variable
x = pde_var(s,[0,1]);
%% Specify the parameters
lam = 9.86;
 if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
 end
%% Define equations
eq_PDE = diff(x,t)==lam*x+diff(x,s,2); % 	PDE: x_{t} = lam*x + x_{ss}
eq_BC = [subs(x,s,0)==0;subs(x,s,1)==0]; %  BCs: x(s=0) = 0,      x(s=1) = 0
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);



% %% Define dependent variables and system variable
% %   PDE: x_{t} = lam*x + x_{ss}                     | lam = 9.86
% %   BCs: x(s=0) = 0,      x(s=1) = 0                |
% x = state('pde');
% pde = sys(); lam = 9.86;
% %% Define equations
% eq_PDE = diff(x,t)==lam*x+diff(x,s,2); % 	PDE: x_{t} = lam*x + x_{ss}
% eq_BC = [subs(x,s,0)==0;subs(x,s,1)==0]; %  BCs: x(s=0) = 0,      x(s=1) = 0
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
