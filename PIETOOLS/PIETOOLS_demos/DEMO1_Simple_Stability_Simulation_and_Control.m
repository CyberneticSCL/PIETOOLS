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
%%
% First, we clear the workspace of any interfering variables
clear all; clc; close all; clear stateNameGenerator


% =============================================
% === Declare the operators of interest

% Declare system parameters
c=1;b=.1;
% Declare independent variables
pvar t s;   
% Declare state, input, and output variables
phi = state('pde',2);   x = state('ode');
w = state('in');        u = state('in');
z = state('out',2);
% Declare the PDE
odepde = sys();
eq_dyn = [diff(x,t,1)==-x+u
          diff(phi,t,1)==[0 1; c 0]*diff(phi,s,1)+[0;s]*w+[0 0;0 -b]*phi];
eq_out= z ==[int([1 0]*phi,s,[0,1]); u];
odepde = addequation(odepde,[eq_dyn;eq_out]);
odepde= setControl(odepde,u);   % set the control signal
bc1 = [0 1]*subs(phi,s,0)==0;   % add the boundary conditions
bc2 = [1 0]*subs(phi,s,1)==x;
odepde = addequation(odepde,[bc1;bc2]);

% Convert to PIE
pie = convert(odepde);


% =============================================
% === Declare the LPI

% Run pre-defined stability executive
[prog, P] = lpiscript(pie,'stability','light');

% Run pre-defined L2-gain executive
[prog, P, gam] = lpiscript(pie,'l2gain','light');

% Run pre-defined controller synthesis executive
[prog, Kval, gam_CL] = lpiscript(pie,'hinf-controller','light');
% Build closed-loop PIE with optimal controller
PIE_CL = closedLoopPIE(pie,Kval);
PIE_CL = pie_struct(PIE_CL);
PIE_CL = initialize(PIE_CL);


% =============================================
% === Simulate the system

% Set PIESIM simulation options
opts.plot = 'no';   % Do not plot the final solution
opts.N = 8;         % Expand using 8 Chebyshev polynomials
opts.tf = 10;       % Simulate up to t = 10;
opts.dt = 0.03;     % Use time step of 10^-2
opts.intScheme = 1; % Time-step using Backward Differentiation Formula (BDF)
ndiff = [0,2,0];    % The PDE state involves 2 first order differentiable state variables   

% Set initial conditions and disturbance for simulation
syms st sx;
uinput.ic.PDE = [0,0]*sx;  
uinput.ic.ODE = 0;  
uinput.u = 0*st;
uinput.w = sin(5*st)*exp(-st); 

% Simulate open-loop PDE and extract solution
[solution,grids] = PIESIM(odepde, opts, uinput, ndiff);
tval = solution.timedep.dtime;
phi1 = reshape(solution.timedep.pde(:,1,:),opts.N+1,[]);
phi2 = reshape(solution.timedep.pde(:,2,:),opts.N+1,[]);
zval = solution.timedep.regulated;
wval = subs(uinput.w,st,tval);

% Simulate closed-loop PIE and extract solution
[solution_CL,grids_CL] = PIESIM(PIE_CL,opts,uinput,ndiff);
tval_CL = solution_CL.timedep.dtime;
phi1_CL = reshape(solution_CL.timedep.pde(:,1,:),opts.N+1,[]);
phi2_CL = reshape(solution_CL.timedep.pde(:,2,:),opts.N+1,[]);
zval_CL =solution_CL.timedep.regulated;
wval_CL = subs(uinput.w,st,tval_CL);


% Plot open-loop state and outputs evolution
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

% Plot closed-loop state and output evolution
figure(3)
surf(tval_CL,grids_CL.phys,phi2_CL,'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
h = colorbar;
colormap jet
box on
ylabel(h,'$|\dot{\mathbf{x}}(t,s)|$','interpreter', 'latex','FontSize',15)
set(gcf, 'Color', 'w');
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$s$','FontSize',15,'Interpreter','latex');
zlabel('$\dot{\mathbf{x}}(t,s)$','FontSize',15,'Interpreter','latex');
title('Closed loop zero-state response with $w(t)=sin(5t)e^{-t}$','Interpreter','latex','FontSize',15);

figure(4);
plot(tval_CL,wval_CL,'k',tval_CL,zval_CL(1,:),'r',tval_CL,zval_CL(2,:),'b','LineWidth',2)
grid on
box on
set(gcf, 'Color', 'w');
legend('$\mathbf{w}(t)$','$\mathbf{r}(t)$','$\mathbf{u}(t)$','Interpreter','latex','FontSize',15)
xlabel('$t$','FontSize',15,'Interpreter','latex');    
ylabel('$\mathbf{r}(t)$','FontSize',15,'Interpreter','latex');
title('Closed loop zero-state response with $w(t)=sin(5t)e^{-t}$','Interpreter','latex','FontSize',15);

% Plot open- and closed-loop output
figure(5);
plot(tval_CL,zval(1,:),'r--',tval_CL,zval_CL(1,:),'r','LineWidth',2)
grid on
box on
set(gcf, 'Color', 'w');
legend('$\mathbf{r}(t)$','$\mathbf{r}_{cl}(t)$','Interpreter','latex','FontSize',15)
xlabel('$t$','FontSize',15,'Interpreter','latex');    
ylabel('$\mathbf{z}(t)$','FontSize',15,'Interpreter','latex');
title('Zero-state response with $w(t)=sin(5t)e^{-t}$','Interpreter','latex','FontSize',15);