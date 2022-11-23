clc; clear;
pvar s t theta; % define independent variables
%% Define dependent variables and system variable
%   PDE: u_{tt} = -c*u_{ssss}                       | c = 0.1
%   BCs: u(s=0) = 0,        u_{ss}(s=1) = 0         |                               Peet 2019 [8] (Example 8.1.0.1)
%        u_{s}(s=0) = 0,    u_{sss}(s=1) = 0        |  
x=state('pde'); 
pde = sys(); 
c = 0.1;
%% Define equations
eq_PDE = diff(x,t,2)==-c*diff(x,s,4); % 	PDE: u_{tt} = -c*u_{ssss}
eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,0)==0;
         subs(diff(x,s,2),s,1)==0;subs(diff(x,s,3),s,1)==0]; %   BCs: u(s=0) = 0, u_{ss}(s=1) = 0, u_{s}(s=0) = 0, u_{sss}(s=1) = 0
%% addequations to pde system; set control inputs/observed inputs, if any
pde = addequation(pde,[eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
