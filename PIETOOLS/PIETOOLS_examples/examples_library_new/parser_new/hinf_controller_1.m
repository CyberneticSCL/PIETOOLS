clc; clear;
%% Define dependent variables and system variable
% % PDE: x_{t} = lam*x + x_{ss} + (s-s^2)u(t) + (s-s^2)w    | lam = 10;
% % BCs: x(s=0) = 0,        x(s=1) = 0
% % Out: z1(t) = int(x(s,t),s,0,1)
% %      z2(t) = u(t)
pvar s t; % define independent variables
pdevar x; outputvar z; inputvar u w; 
lam =10; z.len = 2;
%% Define equations
eq_PDE = [diff(x,t)==lam*x+diff(x,s,2)+(s-s^2)*u+(s-s^2)*w]; %   PDE: x_{t} = lam*x + x_{ss} + (s-s^2)u(t) + (s-s^2)w
eq_BC = [subs(x,s,0)==0;subs(x,s,1)==0]; %   BCs: x(s=0) = 0, x(s=1) = 0
eq_out = [z==[int(x,s,[0,1]);u]]; %   Out: z(t) = [int(x(t,s),s,0,1);u]
%% addequations to pde system; set control inputs/observed inputs, if any
pde = sys('pde',[eq_PDE;eq_BC;eq_out]);
pde = setControl(sys,u);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
