clc; clear;
pvar s t theta; % define independent variables
%% Define dependent variables and system variable
% 	ODE:  xo_{t} = A * xo + Bxr * x_{s}(s=a)        | a = 0    b = 1        
%   PDE:  x_{t}  = lam * x + x_{ss} + Bpv * xo      | lam = pi^2-1          
%   BCs:  x(s=a) = 0,     x(s=b) = 0                | A, Bxr, and Bpv fixed
x=state('ode',4); X = state('pde',2);
pde = sys(); lam = pi^2-1; a = 1;b=1;
A = [-1.2142,  1.9649,  0.2232,  0.5616;
            -1.8042, -0.7260, -0.3479,  5.4355;
            -0.2898,  0.7381, -1.7606,  0.8294;
            -0.9417, -5.3399, -1.0704, -0.7590];
Bpv = [-2.5575 0 1.0368 0;-1.8067 0.4630 1.3621 0];
Bxr = [-1.5368 0;0 0.8871;1.0656 0;1.1882 0];
%% Define equations
eq_ODE=diff(x,t)==A*x+Bxr*diff(subs(X,s,0),s); % ODE: xo_{t} = A * xo + Bxr * x_{s}(s=0)
eq_PDE = diff(X,t)==lam*X+diff(X,s,2)+Bpv*x; % 	PDE: x_{t}  = lam * x + x_{ss} + Bpv * xo
eq_BC = [subs(X,s,0)==0;subs(X,s,1)==0]; %  BCs: x(s=a) = 0,     x(s=b) = 0 
%% addequations to pde system; set control inputs/observed inputs, if any
pde = addequation(pde,[eq_ODE;eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
