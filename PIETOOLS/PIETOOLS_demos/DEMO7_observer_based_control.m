% DEMO7_observer_based_control.m
% See Chapter 11.7 of the manual for a description.
%
% This document illustrates construction of an observer-based state feedback 
% controller. We construct the observer and controller separately
% and couple them to find a closed loop system which is then simulated for 
% different IC and disturbances.
% Specifically, we consider the following system:
% PDE:      \dot{x}(t,s) = (d^2/ds^2) x(t,s) + lam x(t,s) + w(t) + u(t),   s in [0,1];
% Outputs:       z(t)    = [int_{0}^{1} x(t,s) ds; u(t)];
% BCs:                0  = x(t,0) = d/ds x(t,1);
% (unstable for lam > pi^2/4)
% First we convert the above PDE to a PIE of the form:
%   [T \dot{v}](t,s) = [A v](t,s) + [B1 w](t,s) + [B2 u](t,s);
%               z(t) = [C1 v](t)  + [D11 w](t) + [D12 u](t);
%               y(t) = [C2 v](t)  + [D21 w](t) + [D22 u](t)
% We find the observer gain L by solving the LPI
%   min_{gam,P,Z} gam
%   s.t.          P>=0,
%       [-gam*I,           -D11',  -(P*B1+Z*D21)'*T           ] <=0
%       [-D11,             -gam*I, C1                         ]
%       [-T'*(P*B1+Z*D21), C1',    (P*A+Z*C2)'*T+T'*(P*A+Z*C2)]
% where L = P^{-1}*Z.
% Likewise we find a state feedback control u = Kv by solving the LPI
%   min_{gam,P,Z} gam
%   s.t.          P>=0,
%       [ -gam*I            D11      (C1*P+D12*Z)*T'            ]
%       [  D11'            -gam*I    B1'                        ] <= 0
%       [  T*(C1*P+D12*Z)   B1       (A*P+B2*Z)*T'+T*(A*P+B2*Z)']
% where K = Z*P^{-1}. Then, constructing an estimator with gain L, and
% using feedback control u = K*hat{v}, the closed-loop PIE takes the form
%   [T, 0] [\dot{v}(t)   ] = [A,     B2*K  ] [v(t)   ] + [B1   ] w(t)
%   [0, T] [\dot{vhat}(t)]   [-L*C2, A+L*C2] [vhat(t)] + [L*D21]
%
%                [z(t)   ] = [C1, D12*K ] [v   ] + [D11] w(t)
%                [zhat(t)]   [0,  C1    ] [vhat]   [0  ]
% We simulate the open- and closed-loop system response for various
% initial conditions using PIESIM.
%

clc; clear; close all;
echo on

% =============================================
% === Declare the system of interest

% % Declare system as PDE
% Declare independent variables (time and space)
pvar s t
% Declare state, input, and output variables
x = state('pde');       w = state('in');    y = state('out');
z = state('out', 2);    u = state('in');
% Declare the sytem equations
lam = 5;
PDE = sys();
eqs = [diff(x,t) == diff(x,s,2) + lam*x + s*w + s*u;
       z == [int(x,s,[0,1]); u];
       y == subs(x,s,1);
       subs(x,s,0)==0;
       subs(diff(x,s),s,1)==0];
PDE = addequation(PDE,eqs);
PDE = setControl(PDE,u);
PDE = setObserve(PDE,y);
display_PDE(PDE);

% % Convert PDE to PIE
PIE = convert(PDE,'pie');
T = PIE.T;      
A = PIE.A;      B1 = PIE.B1;    B2 = PIE.B2;
C1 = PIE.C1;    D11 = PIE.D11;  D12 = PIE.D12;
C2 = PIE.C2;    D21 = PIE.D21;  D22 = PIE.D22;


% =============================================
% === Declare the LPI

% % Use the predefined Hinf control and estimator functions.
settings = lpisettings('heavy');
[prog_k, Kval, gam_co_val] = lpiscript(PIE,'hinf-controller',settings);
[prog_l, Lval, gam_ob_val] = lpiscript(PIE,'hinf-observer',settings);

% % Construct the closed-loop PIE system
PIE_CL = pie_struct();
PIE_CL.vars = PIE.vars; PIE_CL.dom = PIE.dom;
PIE_CL.T = [T, 0*T; 0*T, T];
PIE_CL.A = [A, B2*Kval; -Lval*C2, A+Lval*C2];   PIE_CL.B1 = [B1; Lval*D21];
PIE_CL.C1 = [C1, D12*Kval; 0*C1, C1];           PIE_CL.D11 = [D11; 0*D11];
PIE_CL = initialize(PIE_CL);

% % Alternatively, uncomment to use predefined functions.
% PIE_CL = closedLoopPIE(PIE,Lval,'observer');
% % Build an operator K_new such that K*vhat = [0,K]*[v; vhat] = K_new*V;
% Kval_new = [0*Kval, Kval];
% % Impose the control law u = K_new*V in the closed-loop observer system
% PIE_CL = closedLoopPIE(PIE_CL,Kval_new,'controller');


% =============================================
% === Simulate the system

% % Declare initial values and disturbance
syms st sx real
uinput1.ic.PDE = [-10*sx; 0];
uinput1.w = 5*exp(-st);
uinput2.ic.PDE = [sin(sx*pi/2); 0];
uinput2.w = 5*sin(pi*st)./(st+eps);

% % Set options for discretization and simulation
opts.plot = 'no';   % don;t plot the final solution
opts.N = 8;         % Expand using 8 Chebyshev polynomials
opts.tf = 10;       % Simulate up to t = 10;
opts.dt = 1e-2;     % Use time step of 10^-2
ndiff = [0,0,1];    % The state involves 1 second order differentiable state variables
ndiff_CL = [0,0,2]; % The closed-loop system involves 2 state variables

% % Perform the actual simulation
% Simulate solution to uncontrolled PIE
[solution_OL_a,grid] = PIESIM(PIE,opts,uinput1,ndiff);
[solution_OL_b,~] = PIESIM(PIE,opts,uinput2,ndiff);
% Simulate solution to controlled PIE
[solution_CL_a,~] = PIESIM(PIE_CL,opts,uinput1,ndiff_CL);
[solution_CL_b,~] = PIESIM(PIE_CL,opts,uinput2,ndiff_CL);
% Extract state and output solution at all times
tval = solution_OL_a.timedep.dtime;
x_OL_a = reshape(solution_OL_a.timedep.pde(:,1,:),opts.N+1,[]);
x_OL_b = reshape(solution_OL_b.timedep.pde(:,1,:),opts.N+1,[]);
x_CL_a = reshape(solution_CL_a.timedep.pde(:,1,:),opts.N+1,[]);
z_CL_a = solution_CL_a.timedep.regulated(1,:);
u_CL_a = solution_CL_a.timedep.regulated(2,:);
x_CL_b = reshape(solution_CL_b.timedep.pde(:,1,:),opts.N+1,[]);
z_CL_b = solution_CL_b.timedep.regulated(1,:);
u_CL_b = solution_CL_b.timedep.regulated(2,:);
w_a = double(subs(uinput1.w,st,tval));
w_b = double(subs(uinput2.w,st,tval));
% Extract state and output estimates at all times
xhat_CL_a = reshape(solution_CL_a.timedep.pde(:,2,:),opts.N+1,[]);
zhat_CL_a = solution_CL_a.timedep.regulated(3,:);
xhat_CL_b = reshape(solution_CL_b.timedep.pde(:,2,:),opts.N+1,[]);
zhat_CL_b = solution_CL_b.timedep.regulated(3,:);


echo off

% % Plot the open-loop state at different times
plot_indcs = 1:5:round(length(tval)/10);
figure('Position',[200 150 1000 450]);
subplot(1,2,1);
surf(tval(plot_indcs),grid.phys,x_OL_a(:,plot_indcs,:),'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
xlabel('$t$','FontSize',15,'Interpreter','latex');  ylabel('$s$','FontSize',15,'Interpreter','latex');  
zlabel('$\mathbf{x}_{ol}(t,s)$','FontSize',15,'Interpreter','latex');
title('$\mathbf{x}(0,s)=-\frac{5}{3}s^3+5s$,\quad $w(t)=e^{-t}$','Interpreter','latex','FontSize',13);
subplot(1,2,2)
surf(tval(plot_indcs),grid.phys,x_OL_b(:,plot_indcs,:),'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
xlabel('$t$','FontSize',15,'Interpreter','latex');  ylabel('$s$','FontSize',15,'Interpreter','latex');  
zlabel('$\mathbf{x}_{ol}(t,s)$','FontSize',15,'Interpreter','latex');
title('$\mathbf{x}(0,s)=-\frac{4}{\pi^2}\sin(\frac{\pi}{2}s)$,\quad $w(t)=\frac{\sin(\pi t)}{t}$','Interpreter','latex','FontSize',13);
sgtitle('Open-loop PDE state','Interpreter','latex','FontSize',15)

% % Plot the closed-loop state at different times
plot_indcs = 1:5:round(length(tval)/10);
figure('Position',[200 150 1000 500]);
subplot(2,2,1);
surf(tval(plot_indcs),grid.phys,x_CL_a(:,plot_indcs,:),'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
xlabel('$t$','FontSize',15,'Interpreter','latex');  ylabel('$s$','FontSize',15,'Interpreter','latex');  
zlabel('$\mathbf{x}_{cl}(t,s)$','FontSize',15,'Interpreter','latex');
title('$\mathbf{x}(0,s)=-\frac{5}{3}s^3+5s$,\quad $w(t)=e^{-t}$','Interpreter','latex','FontSize',13);
subplot(2,2,2)
surf(tval(plot_indcs),grid.phys,x_CL_b(:,plot_indcs,:),'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
xlabel('$t$','FontSize',15,'Interpreter','latex');  ylabel('$s$','FontSize',15,'Interpreter','latex');  
zlabel('$\mathbf{x}_{cl}(t,s)$','FontSize',15,'Interpreter','latex');
title('$\mathbf{x}(0,s)=-\frac{4}{\pi^2}\sin(\frac{\pi}{2}s)$,\quad $w(t)=\frac{\sin(\pi t)}{t}$','Interpreter','latex','FontSize',13);
sgtitle('True and estimated closed-loop PDE state','Interpreter','latex','FontSize',15)
subplot(2,2,3);
surf(tval(plot_indcs),grid.phys,xhat_CL_a(:,plot_indcs,:),'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
xlabel('$t$','FontSize',15,'Interpreter','latex');  ylabel('$s$','FontSize',15,'Interpreter','latex');  
zlabel('$\hat{\mathbf{x}}_{cl}(t,s)$','FontSize',15,'Interpreter','latex');
subplot(2,2,4);
surf(tval(plot_indcs),grid.phys,xhat_CL_b(:,plot_indcs,:),'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
xlabel('$t$','FontSize',15,'Interpreter','latex');  ylabel('$s$','FontSize',15,'Interpreter','latex');  
zlabel('$\hat{\mathbf{x}}_{cl}(t,s)$','FontSize',15,'Interpreter','latex');

% % Plot the closed-loop regulated output versus time
plot_indcs = floor(logspace(0,log(opts.tf/opts.dt)/log(10),50));
tplot = tval(plot_indcs);
colors = {'b','g','m','r','k','c','r','y','o'};     % colors for the plot
figure('Position',[200 150 1000 450]); 
subplot(2,2,1);
hold on
plot(tplot,z_CL_a(1,plot_indcs),[colors{6},'-'],'LineWidth',1.5,'DisplayName','$z_{cl}(t)$');
plot(tplot,zhat_CL_a(1,plot_indcs),[colors{1},'--'],'LineWidth',1.5,'DisplayName','$\hat{z}_{cl}(t)$');
plot(tplot,w_a(plot_indcs),[colors{7},'-.'],'LineWidth',1.0,'DisplayName','$w(t)$');
hold off
subplot(2,2,2);
hold on
plot(tplot,z_CL_b(1,plot_indcs),[colors{6},'-'],'LineWidth',1.5,'DisplayName','$z_{cl}(t)$');
plot(tplot,zhat_CL_b(1,plot_indcs),[colors{1},'--'],'LineWidth',1.5,'DisplayName','$\hat{z}_{cl}(t)$');
plot(tplot,w_b(plot_indcs),[colors{7},'-.'],'LineWidth',1.0,'DisplayName','$w(t)$');
hold off
subplot(2,2,3);
plot(tplot,u_CL_a(1,plot_indcs),[colors{5},'-'],'LineWidth',1.5,'DisplayName','$u_{cl}(t)$');
subplot(2,2,4);
plot(tplot,u_CL_b(1,plot_indcs),[colors{5},'-'],'LineWidth',1.5,'DisplayName','$u_{cl}(t)$');
% Clean up the figure
sgtitle('Closed-loop output $z_{cl}(t)=\int_{0}^{1}\mathbf{x}_{cl}(t,s)ds$, estimate $\hat{z}_{cl}=\int_{0}^{1}\hat{\mathbf{x}}_{cl}(t,s)ds$, and input effort $u_{cl}=\mathcal{K}\hat{\mathbf{x}}_{cl}(t)$','Interpreter','latex','FontSize',15)
ax1 = subplot(2,2,1);   ax1.XScale = 'log';  ax1.XLim = [opts.dt,opts.tf];
legend('Interpreter','latex','Location','northeast','FontSize',10.5);
ylabel('$z_{cl}$','FontSize',15,'Interpreter','latex');
title('$\mathbf{x}(0,s)=-\frac{5}{3}s^3+5s$,\quad $w(t)=e^{-t}$','Interpreter','latex','FontSize',13);
ax2 = subplot(2,2,2);   ax2.XScale = 'log';  ax2.XLim = [opts.dt,opts.tf];
legend('Interpreter','latex','Location','northeast','FontSize',10.5);
ylabel('$z_{cl}$','FontSize',15,'Interpreter','latex');
title('$\mathbf{x}(0,s)=-\frac{4}{\pi^2}\sin(\frac{\pi}{2}s)$,\quad $w(t)=\frac{\sin(\pi t)}{t}$','Interpreter','latex','FontSize',13);
ax3 = subplot(2,2,3);   ax3.XScale = 'log';  ax3.XLim = [opts.dt,opts.tf];
xlabel('$t$','FontSize',15,'Interpreter','latex');  ylabel('$u_{cl}$','FontSize',15,'Interpreter','latex');
ax4 = subplot(2,2,4);   ax4.XScale = 'log';  ax4.XLim = [opts.dt,opts.tf];
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$u_{cl}$','FontSize',15,'Interpreter','latex');
