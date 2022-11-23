% DEMO1_Simple_Stability_Simulation_and_Control_Problem.m
% See Chapter 2 for a full description
%
% This document illustrates, with a simple example, how PIETOOLS can be used to simulate systems
% dynamics, analyse stability and design optimal controllers for the
% system.
%  The Example is a pure transport equation on 1D: 
%  PDE             x_{t} = lam*x + x_{ss} + u(t) +w(t)
%  With BCs        x(s=0) = 0
%                         x(s=1) = 0
% and regulated output z = int_0^1 x ds 

clear all; clc;close all;
pvar t s;syms st st;
x = state('pde'); u = state('in');
w = state('in'); z=state('out',2);
pde = sys();
eq_dyn = diff(x,t)==10*x+diff(x,s,2)+u+w;
eq_out = z==[int(x,s,[0,1]);u];
pde = addequation(pde,[eq_dyn;eq_out]);
eq_bc = [subs(x, s, 0) == 0;subs(x, s, 1) == 0];
pde = addequation(pde,eq_bc);
%% Compute the associated PIE representation, and extract the operators.
PIE = convert(pde,'pie');   
T = PIE.T;
A = PIE.A;      C1 = PIE.C1;    B2 = PIE.B2;
B1 = PIE.B1;    D11 = PIE.D11;  D12 = PIE.D12;
%% Stability Analysis of the system by solving an LPI.
settings = lpisettings('heavy');
[prog, P] = PIETOOLS_stability(PIE,settings);
%% Set options for the discretization and simulation:
opts.plot = 'no';   % Do not plot the final solution
opts.N = 8;         % Expand using 8 Chebyshev polynomials
opts.tf = 10;        % Simulate up to t = 1;
opts.dt = 1e-1;     % Use time step of 10^-3
opts.intScheme=1;   % Time-step using Backward Differentiation Formula (BDF)
ndiff = [0,0,1];    % The PDE state involves 1 second order differentiable state variables
%opts.Norder=1;
%% Simulate the solution to the PDE system without controller.
uinput.ic.PDE = 0;  
uinput.w = 5*exp(-st); % disturbance
[solution,grids] = PIESIM(pde, opts, uinput, ndiff);
%% Extract actual solution at each time step and defining discretized variables.
tval = solution.timedep.dtime;
phi = reshape(solution.timedep.pde(:,1,:),opts.N+1,[]);
zval =solution.timedep.regulated;
wval=subs(uinput.w,st,tval);
%% Plots Open Loop.
figure(1);
surf(tval,grids.phys,phi,'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
h=colorbar ;
colormap jet
box on
ylabel(h,'$|\mathbf{x}(t,s)|$','interpreter', 'latex','FontSize',15)
set(gcf, 'Color', 'w');
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$s$','FontSize',15,'Interpreter','latex');
zlabel('$\mathbf{x}(t,s)$','FontSize',15,'Interpreter','latex');
title('Open loop zero-state response with $w(t)=5e^{-t}$','Interpreter','latex','FontSize',15);

figure(2);
plot(tval,wval,'b',tval,zval,'m','LineWidth',2)
grid on
box on
set(gcf, 'Color', 'w');
legend('$\mathbf{w}(t)$','$\mathbf{z}(t)$','Interpreter','latex','FontSize',15)
xlabel('$t$','FontSize',15,'Interpreter','latex');    
ylabel('$\mathbf{z}(t)$','FontSize',15,'Interpreter','latex');
title('Open loop zero-state response with $w(t)=5e^{-t}$','Interpreter','latex','FontSize',15);

%% Use the predefined Hinf estimator executive function.
settings = lpisettings('heavy',0,'','sedumi');
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
%% Simulate the solution to the PIE with controller
%PIE = convert_PIETOOLS_PDE(odepde);
[solution_CL,grids] = PIESIM(PIE_CL,opts,uinput,ndiff);
%% Extract actual solution at each time step and defining discretized variables.
tval = solution_CL.timedep.dtime;
phi = reshape(solution_CL.timedep.pde(:,1,:),opts.N+1,[]);
%zval =solution_CL.timedep.regulated;
wval=subs(uinput.w,st,tval);

%% Plots Closed Loop.
figure(3);
surf(tval,grids.phys,phi,'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
h=colorbar ;
colormap jet
box on
ylabel(h,'$|\mathbf{x}(t,s)|$','interpreter', 'latex','FontSize',15)
set(gcf, 'Color', 'w');
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$s$','FontSize',15,'Interpreter','latex');
zlabel('$\mathbf{x}(t,s)$','FontSize',15,'Interpreter','latex');
title('Closed loop zero-state response with $w(t)=5e^{-t}$','Interpreter','latex','FontSize',15);

% figure(4);
% plot(tval,wval,'b',tval,zval,'m','LineWidth',2)
% grid on
% box on
% set(gcf, 'Color', 'w');
% legend('$\mathbf{w}(t)$','$\mathbf{z}(t)$','Interpreter','latex','FontSize',15)
% xlabel('$t$','FontSize',15,'Interpreter','latex');    
% ylabel('$\mathbf{z}(t)$','FontSize',15,'Interpreter','latex');
% title('Open loop zero-state response with $w(t)=5e^{-t}$','Interpreter','latex','FontSize',15);