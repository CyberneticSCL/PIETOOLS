% DEMO1_Simple_Stability_Simulation_and_Control_Problem.m
% See Chapter 2 of the manual for a description.
%
% This document illustrates, with a simple example, how PIETOOLS can be 
% used to simulate system dynamics, analyse stability and design optimal 
% controllers for the system.
%  The example is the damped wave equation with dynamic boundary input:
%                xi_{tt} = c*xi_{ss} - b*xi_{t} + s*w(t)
%  using the state phi=(X1,X2) such that X1=xi_{s} and X2=xi_{t}:
%  ODE:     x_{t}  = -x + u(t) 
%  PDEs:    X1_{t} =  X2_{s}
%           X2_{t} = cX1_{s}-bX2_{t}+sw(t)
%  BCs:    X1(s=0) = 0
%          X2(s=1) = x = exp(-t)x(t=0)+ int_0^t exp(-t+tau) u(tau) d tau
%  Output:       z = int_0^1 xi_{s} ds = xi(s=1)-xi(s=0). 
%

clear; clc; close all; clear stateNameGenerator
echo on

% =============================================
% === Declare the system of interest

% % Declare system as PDE
% Declare independent variables
pvar t s;   
% Declare state, input, and output variables
phi = state('pde',2);   x = state('ode');
w = state('in');        u = state('in');
z = state('out',2);
% Declare system equations
c=1;    b=.01;
odepde = sys();
eq_dyn = [diff(x,t,1)==-x+u
          diff(phi,t,1)==[0 1; c 0]*diff(phi,s,1)+[0;s]*w+[0 0;0 -b]*phi];
eq_out= z ==[int([1 0]*phi,s,[0,1]); u];
bc1 = [0 1]*subs(phi,s,0)==0;   % add the boundary conditions
bc2 = [1 0]*subs(phi,s,1)==x;
odepde = addequation(odepde,[eq_dyn;eq_out;bc1;bc2]);
odepde= setControl(odepde,u);   % set the control signal

% % Convert to PIE
pie = convert(odepde);


% =============================================
% === Declare the LPI

% % Run pre-defined stability executive
lpiscript(pie,'stability','light');

% % Run pre-defined L2-gain executive
[~,~,gam] = lpiscript(pie,'l2gain','light');

% % Run pre-defined controller synthesis executive
[prog, Kval, gam_CL] = lpiscript(pie,'hinf-controller','light');
% Build closed-loop PIE with optimal controller
PIE_CL = closedLoopPIE(pie,Kval);


% =============================================
% === Simulate the system

% % Declare initial values and disturbance
syms st sx;
uinput.ic.PDE = [5*sin(2*pi*sx),0];  
uinput.ic.ODE = 0.5;  
uinput.u = 0*st;
uinput.w = sin(5*st)*exp(-st); 

% % Set options for discretization and simulation
opts.plot = 'no';   % don't plot final solution
opts.N = 16;        % expand using 16 Chebyshev polynomials
opts.tf = 9;        % simulate up to t = 9;
opts.dt = 0.03;     % use time step of 3*10^-2
ndiff = [0,2,0];    % PDE state involves 2 first order differentiable state variables   

% % Simulate open-loop PDE and extract solution
[solution,grids] = PIESIM(odepde, opts, uinput, ndiff);
tval = solution.timedep.dtime;
phi1 = reshape(solution.timedep.pde(:,1,:),opts.N+1,[]);
phi2 = reshape(solution.timedep.pde(:,2,:),opts.N+1,[]);
zval = solution.timedep.regulated;
wval = subs(uinput.w,st,tval);

% % Simulate closed-loop PIE and extract solution
[solution_CL,~] = PIESIM(PIE_CL,opts,uinput,ndiff);
tval_CL = solution_CL.timedep.dtime;
phi1_CL = reshape(solution_CL.timedep.pde(:,1,:),opts.N+1,[]);
phi2_CL = reshape(solution_CL.timedep.pde(:,2,:),opts.N+1,[]);
zval_CL = solution_CL.timedep.regulated;
wval_CL = subs(uinput.w,st,tval_CL);


echo off


% % Plot simulated states and regulated outputs against time.
figs_OL = PIESIM_plotsolution(solution,grids,'title','Open-Loop');
figs_CL = PIESIM_plotsolution(solution_CL,grids,'title','Closed-Loop');
% Change titles of plots
fig3 = figs_OL{3};  ax3 = fig3.CurrentAxes;
subtitle(ax3,'Output $z_1(t)$ and Control Effort $u(t)=z_2(t)$','Interpreter','latex','FontSize',13);
fig6 = figs_CL{3};  ax6 = fig6.CurrentAxes;
subtitle(ax6,'Output $z_1(t)$ and Control Effort $u(t)=z_2(t)$','Interpreter','latex','FontSize',13);