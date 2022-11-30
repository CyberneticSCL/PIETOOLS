% DEMO1_Simple_Stability_Simulation_and_Control_Problem.m
% See Chapter 2 of the manual for a description.
%
% This document illustrates, with a simple example, how PIETOOLS can be used to simulate systems
% dynamics, analyse stability and design optimal controllers for the
% system.
%  The example is a pure transport equation on 1D: 
%  PDE                  x_{t}  = lam*x + x_{ss} + u(t) +w(t)
%  With BCs             x(s=0) = 0
%                       x(s=1) = 0
% and regulated output  z      = int_0^1 x ds 
%%
clear all; clc;close all;
pvar t s;syms st sx;
%% Declaring the System
% define the parameters
c=1;b=.1;
%%%%%%%%%%%%%%%%%%%%%%
% create the equation variables
phi=state('pde',2);x=state('ode');
w= state('in');u=state('in');
z=state('out',2);
% create the system
odepde = sys();
% add the dynamic equations
eq_dyn = [diff(x,t,1) == -x+u
                diff(phi,t,1)==[0 1; c 0]*diff(phi,s,1)+[0;s]*w+[0 0;0 -b]*phi];
eq_out= z ==[int([1 0]*phi,s,[0,1])
                                                    u];
odepde = addequation(odepde,[eq_dyn;eq_out]);
% add the boundary conditions
bc1 = [0 1]*subs(phi,s,0) == 0;
bc2 = [1 0]*subs(phi,s,1) == x;
odepde = addequation(odepde,[bc1;bc2]);
% set the control signal
odepde= setControl(odepde,[u]);
%% Set options for the discretization and simulation:
opts.plot = 'no';   % Do not plot the final solution
opts.N = 8;         % Expand using 8 Chebyshev polynomials
opts.tf = 10;        % Simulate up to t = 1;
opts.dt = 1e-2;     % Use time step of 10^-3
opts.intScheme=1;   % Time-step using Backward Differentiation Formula (BDF)
ndiff = [0,2,0];    % The PDE state involves 2 first order differentiable state variables   
%opts.Norder=1;
%% Simulate the solution to the PDE system without controller.
uinput.ic.PDE = [0,0];  
uinput.ic.ODE = 0;  
uinput.u=0;
uinput.w = sin(5*st)*exp(-st); % disturbance
[solution,grids] = PIESIM(odepde, opts, uinput, ndiff);
%% Extract actual solution at each time step and defining discretized variables.
tval = solution.timedep.dtime;
phi1 = reshape(solution.timedep.pde(:,1,:),opts.N+1,[]);
phi2 = reshape(solution.timedep.pde(:,2,:),opts.N+1,[]);
zval =solution.timedep.regulated;
wval=subs(uinput.w,st,tval);
%% Plots of open-loop system.
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
legend('$\mathbf{w}(t)$','$\mathbf{r}(t)$','$\mathbf{u}(t)$','Interpreter','latex','FontSize',15)
xlabel('$t$','FontSize',15,'Interpreter','latex');    
ylabel('$\mathbf{r}(t)$','FontSize',15,'Interpreter','latex');
title('Open loop zero-state response with $w(t)=sin(5t)5e^{-t}$','Interpreter','latex','FontSize',15);
%% Stability Analysis of the system by solving an LPI.
% compute the associated PIE representation, and extract the operators.
PIE = convert(odepde,'pie');   
T = PIE.T;
A = PIE.A;      C1 = PIE.C1;    B2 = PIE.B2;
B1 = PIE.B1;    D11 = PIE.D11;  D12 = PIE.D12;
% call the executive with chosen settings
settings = lpisettings('heavy');
[prog, P] = PIETOOLS_stability(PIE,settings);
%% Hinf gain of the open-loop system.
[prog, P, gamma] = PIETOOLS_Hinf_gain(PIE,settings);

%% Use the predefined Hinf estimator executive function.
[prog, Kval, gam_val] = PIETOOLS_Hinf_control(PIE, settings);

%% Closed loop with the sinthesized controller.
    T_CL = T;
    A_CL = A+B2*Kval;   B_CL = B1;
    C_CL = C1+D12*Kval; D_CL = D11;
%% Declare the PIE.
    PIE_CL = pie_struct();
    PIE_CL.vars = PIE.vars;
    PIE_CL.dom = PIE.dom;
    PIE_CL.T = T_CL;
    PIE_CL.A = A_CL;        PIE_CL.B1 = B_CL;
    PIE_CL.C1 = C_CL;       PIE_CL.D11 = D_CL;
    PIE_CL = initialize(PIE_CL);
%     %% Stability Analysis of the system by solving an LPI.
% settings = lpisettings('heavy');
% [prog, P] = PIETOOLS_stability(PIE_CL,settings);
%% Simulate the solution to the PIE with controller
%PIE = convert_PIETOOLS_PDE(odepde);
[solution_CL,grids] = PIESIM(PIE_CL,opts,uinput,ndiff);
%% Extract actual solution at each time step and defining discretized variables.
tval = solution_CL.timedep.dtime;
phi1 = reshape(solution_CL.timedep.pde(:,1,:),opts.N+1,[]);
phi2 = reshape(solution_CL.timedep.pde(:,2,:),opts.N+1,[]);
zval_cl =solution_CL.timedep.regulated;
wval=subs(uinput.w,st,tval);

%% Plots Closed Loop.
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
%%
figure(5);
plot(tval,zval(1,:),'r',tval,zval_cl(1,:),'b','LineWidth',2)
grid on
box on
set(gcf, 'Color', 'w');
legend('$\mathbf{r}(t)$','$\mathbf{r}_{cl}(t)$','Interpreter','latex','FontSize',15)
xlabel('$t$','FontSize',15,'Interpreter','latex');    
ylabel('$\mathbf{z}(t)$','FontSize',15,'Interpreter','latex');
title('Closed loop zero-state response with $w(t)=sin(5t)e^{-t}$','Interpreter','latex','FontSize',15);