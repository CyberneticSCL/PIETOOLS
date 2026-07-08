% Timoshenko Beam
%   PDE: r*aa * w_{tt} = k*aa*g * (-phi_{s} + w_{ss})                               Peet 2019 [8] (Example 8.1.0.2)
%        r*II * phi_{tt} = E*II * phi_{ss}  + k*aa*g * (w_{s} - phi)-c
%        phi_{t}  
%   BCs: phi(s=0) = 0,      phi_{s}(s=1) = 0    
%        w(s=0) = 0,        w_{s}(s=1) - phi(s=1) = 0           
% 
%   Use states:                                     | k = 1                 Hyperbolic implementation | (PDE stable without damping c=0)
%        x1 = w_{t},   x2 = k*aa*g * (w_{s}-phi),   | aa = 1 
%        x3 = phi_{t}-c phi_{t}, x4 = E*II * phi_{s}.         | II = 1   
%       =>                                          | g = 1
%   PDE: x1_{t} = (1/r/aa) * x2_{s}                 | E = 1   (PDE stable and Finite-energy PDE stable with damping c >0)
%        x2_{t} = k*aa*g * x1_{s} - k*aa*g * x3     | r = 1
%        x3_{t} = (1/r/II)*x2 + (1/r/II)*x4_{s}     | 
%        x4_{t} = E*II * x3_{s}                     |
%   BCs: x1(0) = 0,         x2(1) = 0               |
%        x3(0) = 0,         x4(1) = 0               |
%                                                   |
clc; clear; clear stateNameGenerator
pvar s t; % define independent variables
%% Specify the parameters 
k = 1;    aa = 1;    II = 1;    g = 1;    c=0;
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%% Define dependent variables and system variable    
pde_var x1 x2 x3 x4
x1.vars = s;    x2.vars = s;    x3.vars = s;    x4.vars = s;
%% Define equations
%   PDE: x1_{t} = (1/aa) * x2_{s}                
%        x2_{t} = k*aa*g * x1_{s} - k*aa*g * x3   
%        x3_{t} = (1/II)*x2 + (1/II)*x4_{s}  -c x_3   
%        x4_{t} = II * x3_{s}
eq_PDE = [diff(x1,t)==(1/aa)*diff(x2,s); 
          diff(x2,t)==k*aa*g*diff(x1,s)-k*aa*g*x3;
          diff(x3,t)==(1/II)*x2+(1/II)*diff(x4,s)-c*x3;
          diff(x4,t)==II*diff(x3,s)];
eq_BC = [subs(x1,s,0)==0;subs(x3,s,0)==0;
         subs(x2,s,1)==0;subs(x4,s,1)==0]; %   BCs:BCs: x1(0) = 0, x2(1) = 0, x3(0) = 0, x4(1) = 0 
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);



% %% Define dependent variables and system variable
% %   PDE: x1_{t} = (1/aa) * x2_{s}                
% %        x2_{t} = k*aa*g * x1_{s} - k*aa*g * x3   
% %        x3_{t} = (1/II)*x2 + (1/II)*x4_{s}     
% %        x4_{t} = II * x3_{s}                   
% %   BCs: x1(0) = 0,         x2(1) = 0           
% %        x3(0) = 0,         x4(1) = 0           
% x1=state('pde');x2=state('pde');x3=state('pde');x4=state('pde'); 
% pde = sys(); 
% k = 1;    aa = 1;    II = 1;    g = 1;    
% %% Define equations
% %   PDE: x1_{t} = (1/aa) * x2_{s}                
% %        x2_{t} = k*aa*g * x1_{s} - k*aa*g * x3   
% %        x3_{t} = (1/II)*x2 + (1/II)*x4_{s}     
% %        x4_{t} = II * x3_{s}
% eq_PDE = [diff(x1,t)==(1/aa)*diff(x2,s); 
%           diff(x2,t)==k*aa*g*diff(x1,s)-k*aa*g*x3;
%           diff(x3,t)==(1/II)*x2+(1/II)*diff(x4,s);
%           diff(x4,t)==II*diff(x3,s)];
% eq_BC = [subs(x1,s,0)==0;subs(x3,s,0)==0;
%          subs(x2,s,1)==0;subs(x4,s,1)==0]; %   BCs:BCs: x1(0) = 0, x2(1) = 0, x3(0) = 0, x4(1) = 0 
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
