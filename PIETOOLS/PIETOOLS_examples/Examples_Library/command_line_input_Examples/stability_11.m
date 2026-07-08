% 	PDE:  x_{t} = Cm*x + (1/R)*x_{ss}               | R = 2.9               (stable for R<2.7 (not tight))
%   BCs:  x(s=0) = 0,     x(s=1) = 0                | Cm = [1, 1.5; 5, 0.2]         Ahmadi 2014 [6] Example D
%                                                                         (Exponentially PDE stable )
%                                                                         (Finite Energy PDE stable ) 
clc; clear; clear stateNameGenerator
pvar s t; % define independent variables
% Define dependent variables and system variable    
x=pde_var(2,s,[0,1]); 
%% Specify the parameters
R = 2,9;Cm = [1, 1.5; 5, 0.2];
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%% Define equations
eq_PDE = diff(x,t)==Cm*x+(1/R)*diff(x,s,2); % 	PDE: x_{t} = Cm*x + (1/R)*x_{ss}
eq_BC = [subs(x,s,0)==0;subs(x,s,1)==0]; %  BCs: x(s=0) = 0,     x(s=1) = 0 
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);



% %% Define dependent variables and system variable
% % 	PDE:  x_{t} = Cm*x + (1/R)*x_{ss}               | R = 2.7
% %   BCs:  x(s=0) = 0,     x_{s}(s=1) = 0            | Cm = [1, 1.5; 5, 0.2]     
% x=state('pde',2); 
% pde = sys(); 
% R = 2.7;Cm = [1, 1.5; 5, 0.2];
% %% Define equations
% eq_PDE = diff(x,t)==Cm*x+(1/R)*diff(x,s,2); % 	PDE: x_{t} = Cm*x + (1/R)*x_{ss}
% eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,1)==0]; %  BCs: x(s=0) = 0,     x_{s}(s=1) = 0 
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
