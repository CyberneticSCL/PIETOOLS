%% Introduction
% DEMO1_Simple_Stability_Simulation_and_Control_Problem.m
% See Chapter 2 of the manual for a description.
%
% This document illustrates, with a simple example, how PIETOOLS can be used to simulate systems
% dynamics, analyse stability and design optimal controllers for the
% system.
%  The example is the damped wave equation with dynamic boundary input:
%                xi_{tt}=c xi_{ss}-b xi_{t}+sw(t)
%  using the state phi=(X1,X2) such that X1=xi_{s} and X2=xi_{t}:
%  ODE                 x_{t}  = -x + u(t) 
%  PDEs:              X1_{t} =  X2_{s}
%                          X2_{t} = cX1_{s}-bX2_{t}+sw(t)
%  With BCs         X1(s=0) = 0
%                          X2(s=1) = x = exp(-t)x(t=0)+ int_0^t exp(-t+tau) u(tau) d tau
%                  
% and regulated output  z= int_0^1 xi_{s} ds = xi(s=1)-xi(s=0). 
%%
% First, we clear the workspace of any interfering variables
clear all; clc;close all;
%% Declaring variables for symbolic manipulation
% Define the parameters and variables: Here we declare symbolic variables
% t (st) and s (sx) that stand for time and space, respectively, in
% polynomial (symbolic) object format. Then, we define states {phi, x},
% inputs {w, u}, and output {z}, which will be used to define the PDE
% system described earlier.
c=1;b=.1;
pvar t s; syms st sx;
phi=state('pde',2);x=state('ode');
w= state('in');u=state('in');
z=state('out',2);
%% Create the system 
% Here we use a sys() class object to store the equations that describe
% the PDE model. First, we initialize the object and then equations are
% added to the object using addequation() function. Finally, setControl()
% function is used designate u as a control input.
odepde = sys();
eq_dyn = [diff(x,t,1)==-x+u
                diff(phi,t,1)==[0 1; c 0]*diff(phi,s,1)+[0;s]*w+[0 0;0 -b]*phi];
eq_out= z ==[int([1 0]*phi,s,[0,1])
                                 u];
odepde = addequation(odepde,[eq_dyn;eq_out]);
odepde= setControl(odepde,[u]); % set the control signal
bc1 = [0 1]*subs(phi,s,0)==0; % add the boundary conditions
bc2 = [1 0]*subs(phi,s,1)==x;
odepde = addequation(odepde,[bc1;bc2]);
%% Simulating the Open Loop system
% Here we simulate the PDE system without any control input to see the
% open-loop behaviour of the system using PIESIM. The following options are
% used to define simulation parameters, initial conditions and disturbance
% inputs for the simulation. 
opts.plot = 'no';   % Do not plot the final solution
opts.N = 8;         % Expand using 8 Chebyshev polynomials
opts.tf = 10;        % Simulate up to t = 1;
opts.dt = 1e-2;     % Use time step of 10^-3
opts.intScheme=1;   % Time-step using Backward Differentiation Formula (BDF)
ndiff = [0,2,0];    % The PDE state involves 2 first order differentiable state variables   

% To simulate the solution to the PDE system without controller we can use
% the following code. 
uinput.ic.PDE = [0,0]*sx;  
uinput.ic.ODE = 0;  
uinput.u=0*st;
uinput.w = sin(5*st)*exp(-st); 
[solution,grids] = PIESIM(odepde, opts, uinput, ndiff);

% Extract actual solution at each time step and defining discretized variables.
tval = solution.timedep.dtime;
phi1 = reshape(solution.timedep.pde(:,1,:),opts.N+1,[]);
phi2 = reshape(solution.timedep.pde(:,2,:),opts.N+1,[]);
zval =solution.timedep.regulated;
wval=subs(uinput.w,st,tval);

% Plots of open-loop system.
figure(1);
surf(tval,grids.phys,phi2,'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
h=colorbar ;
colormap jet
box on
ylabel(h,'$|\dot{\mathbf{x}}(t,s)|$','interpreter', 'latex','FontSize',15)
set(gcf, 'Color', 'w');
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$s$','FontSize',15,'Interpreter','latex');
zlabel('$\dot{\mathbf{x}}(t,s)$','FontSize',15,'Interpreter','latex');
title('Open loop zero-state response with $w(t)=sin(5t)e^{-t}$','Interpreter','latex','FontSize',15);

figure(2);
plot(tval,wval,'k',tval,zval(1,:),'r','LineWidth',2)
grid on
box on
set(gcf, 'Color', 'w');
legend('$\mathbf{w}(t)$','$\mathbf{r}(t)$','Interpreter','latex','FontSize',15)
xlabel('$t$','FontSize',15,'Interpreter','latex');    
ylabel('$\mathbf{r}(t)$','FontSize',15,'Interpreter','latex');
title('Open loop zero-state response with $w(t)=sin(5t)e^{-t}$','Interpreter','latex','FontSize',15);
%% Stability Analysis of the system
% To use LPIs for the test of internal stability, we have to first convert
% the PDE representation to a PIE representation. This is performed using
% the following code. 
PIE = convert(odepde,'pie'); 
%PIE = PIE.params;   
T = PIE.T;
A = PIE.A;      C1 = PIE.C1;    B2 = PIE.B2;
B1 = PIE.B1;    D11 = PIE.D11;  D12 = PIE.D12;

% Pick parameters for the LPI optimization problem and perform stability
% test.
settings = lpisettings('light');
[prog, P] = lpiscript(PIE,'stability',settings);

% Other LPI tests:
% Similarly, we can search for an input-to-output L2-gain bound using LPIs
% as shown below.
[prog, P, gamma] = lpiscript(PIE,'l2gain',settings);

%% Find Hinf-optimal controller
% Next, we wish to find a controller that improves the input-to-output
% performance. For that, we can solve the Hinf-optimal controller LPI using
% the following code. The function, if the optimization problem is
% successfully solved, returns the controller, Kval and the L2-gain metric
% for the closed loop system, gam_val.
[prog, Kval, gam_val] = lpiscript(PIE,'hinf-controller',settings);

%% Constructing closed-loop system and simulation
% Now, we construct the closed-loop system using the controller obtained by
% solving the above LPI and then re-runs the simulations to see if there
% was any improvement in the performance.
PIE_CL = closedLoopPIE(PIE,Kval);
PIE_CL = pie_struct(PIE_CL);
PIE_CL = initialize(PIE_CL);

% Simulate the solution to the PIE with controller
[solution_CL,grids] = PIESIM(PIE_CL,opts,uinput,ndiff);

%% Plotting the closed loop response
% Having found an optimal controller, constructed closed loop system, and
% running PIESIM to simulate, we can now extract the solution and plot it
% against the open-loop response to see the benefit in using the
% controller. Firstly, we notice that the L2-gain bound has significantly
% lowered (from 5.2 to 0.78). Then, looking at the output in presence of a
% disturbance, we see that the effect of disturbance on this neutrally
% stable system is reduced.
tval = solution_CL.timedep.dtime;
phi1 = reshape(solution_CL.timedep.pde(:,1,:),opts.N+1,[]);
phi2 = reshape(solution_CL.timedep.pde(:,2,:),opts.N+1,[]);
zval_cl =solution_CL.timedep.regulated;
wval=subs(uinput.w,st,tval);

% Plots Closed Loop.
figure(3)
surf(tval,grids.phys,phi2,'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
h=colorbar ;
colormap jet
box on
ylabel(h,'$|\dot{\mathbf{x}}(t,s)|$','interpreter', 'latex','FontSize',15)
set(gcf, 'Color', 'w');
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$s$','FontSize',15,'Interpreter','latex');
zlabel('$\dot{\mathbf{x}}(t,s)$','FontSize',15,'Interpreter','latex');
title('Closed loop zero-state response with $w(t)=sin(5t)e^{-t}$','Interpreter','latex','FontSize',15);

figure(4);
plot(tval,wval,'k',tval,zval_cl(1,:),'r',tval,zval_cl(2,:),'b','LineWidth',2)
grid on
box on
set(gcf, 'Color', 'w');
legend('$\mathbf{w}(t)$','$\mathbf{r}(t)$','$\mathbf{u}(t)$','Interpreter','latex','FontSize',15)
xlabel('$t$','FontSize',15,'Interpreter','latex');    
ylabel('$\mathbf{r}(t)$','FontSize',15,'Interpreter','latex');
title('Closed loop zero-state response with $w(t)=sin(5t)e^{-t}$','Interpreter','latex','FontSize',15);

figure(5);
plot(tval,zval(1,:),'r--',tval,zval_cl(1,:),'r','LineWidth',2)
grid on
box on
set(gcf, 'Color', 'w');
legend('$\mathbf{r}(t)$','$\mathbf{r}_{cl}(t)$','Interpreter','latex','FontSize',15)
xlabel('$t$','FontSize',15,'Interpreter','latex');    
ylabel('$\mathbf{z}(t)$','FontSize',15,'Interpreter','latex');
title('Zero-state response with $w(t)=sin(5t)e^{-t}$','Interpreter','latex','FontSize',15);