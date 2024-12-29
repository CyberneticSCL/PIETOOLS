
% This document illustrates how PIETOOLS can be used to simulate system
% dynamics, analyse stability and design optimal controllers.
% We refer to Chapter 2 of the manual for more context on the codes.
%  
% For this example, we consider the damped wave equation with dynamic 
% boundary input,
%                xi_{tt}=c xi_{ss}-b xi_{t}+sw(t),
% and with regulated output
%               r(t) = int_0^1 xi_{s} ds = xi(s=1)-xi(s=0).
% Introducing the state phi=(phi1,phi2) with phi1=xi_{s} and phi2=xi_{t},
% this wave equation can be expressed as an ODE-PDE system
%  ODE                 x_{t} = -x + u(t) 
%  PDEs:            phi1_{t} = phi2_{s}
%                   phi2_{t} = c phi1_{s} - b phi2_{t}+sw(t)
%  With BCs        phi1(s=0) = 0
%                  phi2(s=1) = x = exp(-t)x(t=0)+ int_{0}^{t} exp(-t+tau) u(tau) d tau
%  and outputs         z1(t) = int_{0}^{1} phi1(t,s) ds
%                      z2(t) = u(t)
% where now we add an output z2(t)=u(t) so that, in constructing an optimal
% control, we also minimize the control effort.
% DB, 12/29/2024: Use pde_var objects instead of sys and state


%% 2.2.1 Declare the ODE-PDE model
% Declare an ODE state variable x, a vector-valued 1D PDE state variable
% phi with 2 elements, and a finite-dimensional exogenous input w,
% controlled input u, and regulated output z, with z containing 2 elements. 
clc; clear; close all; clear stateNameGenerator;
% Declare an ODE-PDE system structure 'sys' describing the system
% Declare independent variables
pvar t s;   
% Declare state, input, and output variables
phi = pde_var('state',2,s,[0,1]);   x = pde_var('state',1,[],[]);
w = pde_var('input',1);     u=pde_var('control',1);   
z = pde_var('output',2);
% Declare system parameters
c=1;    b=.01;
% Declare ODE-PDE equations
%   d/dt      x(t) = -x(t) + u(t);
%   d/dt phi1(t,s) = phi2_{ss}(t,s);
%   d/dt phi2(t,s) = c*phi1_{ss}(t,s) + s*w(t) -b*phi2(t,s);
eq_dyn = [diff(x,t,1)==-x+u
          diff(phi,t,1)==[0 1; c 0]*diff(phi,s,1)+[0;s]*w+[0 0;0 -b]*phi];
% Declare output equations
%            z1(t) = int_{0}^{1} phi1(t,s)ds;
%            z2(t) = u(t);
eq_out= z==[int([1 0]*phi,s,[0,1]); 
            u];

%% 2.2.2 Declare the boundary conditions
bc1 = [0 1]*subs(phi,s,0) == 0;         % phi2(t,0) = 0;
bc2 = [1 0]*subs(phi,s,1) == x;         % phi1(t,1) = x(t);
% create the PDE system
odepde = [eq_dyn;eq_out;bc1;bc2];

%% 2.2.3: Simulate the ODE-PDE system

% Declare simulation settings.
opts.plot = 'no';   % don't plot final solution
opts.N = 8;         % expand solution using 8 Chebyshev polynomials
opts.tf = 15;       % simulate up to t = 15
opts.dt = 3*1e-2;   % use time step of 3*10^-2
opts.intScheme = 1; % use backward differentiation formula for time integration

% Declare initial values and values of disturbance and input.
syms st sx real                 % symbolic time and space variables
uinput.ic.PDE = [0,0];          % initial values of both PDE state variables
uinput.ic.ODE = 0;              % initial value of ODE state variable
uinput.u = 0;                   % value of control as function of time
uinput.w = sin(5*st)*exp(-st);  % value of disturbance as function of time

% Declare number of PDE state variables differentiable up to each order:
%   0 0th-order, 2 1st-order, and 0 2nd-order differentiable states.
ndiff = [0,2,0]; 

% Simulate the ODE-PDE system using specified settings.
[solution,grids] = PIESIM(odepde, opts, uinput, ndiff);

% Extract values of the second state variables and thet output at each time step.
tval = solution.timedep.dtime;
phi2 = reshape(solution.timedep.pde(:,2,:),opts.N+1,[]);
zval = solution.timedep.regulated;
wval = subs(uinput.w,st,tval);

% Plot evolution of the PDE states.
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

% Plot evolution of the disturbance and regulated output.
figure(2);
plot(tval,wval,'k',tval,zval(1,:),'r','LineWidth',2)
grid on
box on
set(gcf, 'Color', 'w');
legend('$\mathbf{w}(t)$','$\mathbf{r}(t)$','Interpreter','latex','FontSize',15)
xlabel('$t$','FontSize',15,'Interpreter','latex');    
ylabel('$\mathbf{r}(t)$','FontSize',15,'Interpreter','latex');
title('Open loop zero-state response with $w(t)=sin(5t)e^{-t}$','Interpreter','latex','FontSize',15);


%% 2.2.4: Analyse stability and construct a controller

% Convert the ODE-PDE to a PIE.
PIE = convert(odepde,'pie');

% Verify stability of the system by runnig the stability executive.
settings = lpisettings('heavy');    % use heavy settings: more expensive but less conservative
lpiscript(PIE,'stability',settings);

% Compute an upper bound gam on the H_infty norm or L2-gain of the system,
% so that ||z(.)||_{L2} <= gam*||w(.)||_{L2}.
[~,~, gam] = lpiscript(PIE,'l2gain',settings);

% Compute a feedback gain operator K that minimizes the H_infty norm of the
% closed-loop system with state feedback u = K*v, for fundamental (PIE)
% state v.
[~, Kval, gam_CL] = lpiscript(PIE,'hinf-controller','light');

% Now, we construct the closed-loop system using the controller obtained by
% solving the above LPI and then re-runs the simulations to see if there
% was any improvement in the performance.
PIE_CL = closedLoopPIE(PIE,Kval);
PIE_CL = pie_struct(PIE_CL);
PIE_CL = initialize(PIE_CL);


%% -.-.- Simulate the closed-loop ODE-PDE system
% Simulate the solution to the PIE with controller
[solution_CL,grids_CL] = PIESIM(PIE_CL,opts,uinput,ndiff);

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
zval_cl = solution_CL.timedep.regulated;
wval = subs(uinput.w,st,tval);

% Plot closed-loop evolution of second state variable.
figure(3)
surf(tval,grids_CL.phys,phi2,'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
h=colorbar ;
colormap jet
box on
ylabel(h,'$|\dot{\mathbf{x}}(t,s)|$','interpreter', 'latex','FontSize',15)
set(gcf, 'Color', 'w');
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$s$','FontSize',15,'Interpreter','latex');
zlabel('$\dot{\mathbf{x}}(t,s)$','FontSize',15,'Interpreter','latex');
title('Closed loop zero-state response with $w(t)=sin(5t)e^{-t}$','Interpreter','latex','FontSize',15);

% Plot closed-loop evolution of disturbance w and regulated output r(t).
figure(4);
plot(tval,wval,'k',tval,zval_cl(1,:),'r',tval,zval_cl(2,:),'b','LineWidth',2)
grid on
box on
set(gcf, 'Color', 'w');
legend('$\mathbf{w}(t)$','$\mathbf{r}(t)$','$\mathbf{u}(t)$','Interpreter','latex','FontSize',15)
xlabel('$t$','FontSize',15,'Interpreter','latex');    
ylabel('$\mathbf{r}(t)$','FontSize',15,'Interpreter','latex');
title('Closed loop zero-state response with $w(t)=sin(5t)e^{-t}$','Interpreter','latex','FontSize',15);

% Compare evolution of open-loop and closed-loop regulated output r(t).
figure(5);
plot(tval,zval(1,:),'r--',tval,zval_cl(1,:),'r','LineWidth',2)
grid on
box on
set(gcf, 'Color', 'w');
legend('$\mathbf{r}(t)$','$\mathbf{r}_{cl}(t)$','Interpreter','latex','FontSize',15)
xlabel('$t$','FontSize',15,'Interpreter','latex');    
ylabel('$\mathbf{z}(t)$','FontSize',15,'Interpreter','latex');
title('Zero-state response with $w(t)=sin(5t)e^{-t}$','Interpreter','latex','FontSize',15);