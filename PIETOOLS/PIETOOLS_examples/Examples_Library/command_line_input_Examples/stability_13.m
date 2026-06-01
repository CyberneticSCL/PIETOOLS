clc; clear; clear stateNameGenerator
pvar s t; % define independent variables
%% Define dependent variables and system variable
% % This case uses PDE states x1 = u_{t}, x2 = u_{ss}, so that the system becomes
% % PDE        x1_{t} = -c * x2_{ss}
% %            x2_{t} = x1_{ss}
% % with BCs   x1(0) = 0, x2(1) = 0, x1_{s}(0) = 0, x2_{s}(1) = 0
% Construct the PDE.
x1=pde_var(s,[0,1]); x2=pde_var(s,[0,1]); 
c = 0.1;
%% Define equations
eq_PDE = [diff(x1,t,1)==-c*diff(x2,s,2); % PDE        x1_{t} = -c * x2_{ss}
                   diff(x2,t,1)==diff(x1,s,2);]%           x2_{t} = x1_{ss}
eq_BC = [subs(x1,s,0)==0;subs(diff(x1,s,1),s,0)==0;% with BCs   x1(0) = 0, x1_{s}(0) = 0, 
         subs(x2,s,1)==0;subs(diff(x2,s,1),s,1)==0];  % x2(1) = 0, x2_{s}(1) = 0
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);



% %% Define dependent variables and system variable
% %   PDE: u_{tt} = -c*u_{ssss}                       | c = 0.1
% %   BCs: u(s=0) = 0,        u_{ss}(s=1) = 0         |                               Peet 2019 [8] (Example 8.1.0.1)
% %        u_{s}(s=0) = 0,    u_{sss}(s=1) = 0        |  
% x=state('pde'); 
% pde = sys(); 
% c = 0.1;
% %% Define equations
% eq_PDE = diff(x,t,2)==-c*diff(x,s,4); % 	PDE: u_{tt} = -c*u_{ssss}
% eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,0)==0;
%          subs(diff(x,s,2),s,1)==0;subs(diff(x,s,3),s,1)==0]; %   BCs: u(s=0) = 0, u_{ss}(s=1) = 0, u_{s}(s=0) = 0, u_{sss}(s=1) = 0
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
