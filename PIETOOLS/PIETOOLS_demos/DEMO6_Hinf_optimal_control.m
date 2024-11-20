% DEMO6_Hinf_optimal_control.m
% See Chapter 11.6 of the manual for a description.
%
% This document illustrates how an Hinfty optimal controller can be designed
% for a PDE using the PIE/LPI framework.
% Specifically, we consider the following system:
% PDE:      \dot{x}(t,s) = (d^2/ds^2) x(t,s) + lam x(t,s) + w(t) + u(t),   s in [0,1];
% Outputs:       z(t)    = [int_{0}^{1} x(t,s) ds + w(t); u(t)];
% BCs:                0  = x(t,0) = d/ds x(t,1);
% (unstable for lam > pi^2/4)
% Letting v:=(d^2/ds^2)x, We derive an equivalent PIE of the form:
%   [T \dot{v}](t,s) = [A v](t,s) + [B1 w](t,s) + [B2 u](t,s);
%               z(t) = [C1 v](t)  + [D11 w](t) + [D12 u](t);
% Using a state feedback control u = K*v we get the closed loop PIE
%   [T \dot{v}](t,s) = [(A + B2*K) v](t,s) + [B1 w](t,s);
%               z(t) = [(C1 + D12*K) v](t) + [D11 w](t)
% To compute an operator K that minimizes the L2 gain from disturbance w to
% the output z, we solve the LPI
%   min_{gam,P,Z}   gam
%   s.t.            P>=0,
%       [ -gam*I            D11      (C1*P+D12*Z)*T'            ]
%       [  D11'            -gam*I    B1'                        ] <= 0
%       [  T*(C1*P+D12*Z)   B1       (A*P+B2*Z)*T'+T*(A*P+B2*Z)']
%
% Then, using K = Z*P^{-1}, the closed-loop system with u = K*v satisfies
%   ||z||_{L2}/||w||_{L2} <= gam
% We manually declare this LPI here, but it can also be solved using the
% "PIETOOLS_Hinf_control" executive file.
% We simulate the open loop and closed loop response of the PDE for various
% IC using PIESIM.
%

clc; clear; close all;
echo on

% =============================================
% === Declare the system of interest

% % Declare system as PDE
% Declare independent variables (time and space)
pvar s t
% Declare state, input, and output variables
x = state('pde');       w = state('in');
z = state('out', 2);    u = state('in');
% Declare the sytems equations
lam = 5;
PDE = sys();
eqs = [diff(x,t) == diff(x,s,2) + lam*x + s*w + s*u;
       z == [int(x,s,[0,1])+w; u];
       subs(x,s,0)==0;
       subs(diff(x,s),s,1)==0];
PDE = addequation(PDE,eqs);
% Set u as constrolled input
PDE = setControl(PDE,u);
display_PDE(PDE);

% % Convert PDE to PIE
PIE = convert(PDE,'pie');
T = PIE.T;
A = PIE.A;      C1 = PIE.C1;    B2 = PIE.B2;
B1 = PIE.B1;    D11 = PIE.D11;  D12 = PIE.D12;


% =============================================
% === Declare the LPI

% % Initialize LPI program
prog = lpiprogram(PIE.vars,PIE.dom);

% % Declare decision variables:
% %   gam \in \R,     P:L2-->L2,    Z:\R-->L2
% Scalar decision variable
dpvar gam
prog = lpidecvar(prog,gam);
% Positive operator variable P>=0
Pdim = T.dim(:,1);
Pdeg = {2,[1,1,2],[1,1,2]};
opts.sep = 0;
[prog,P] = poslpivar(prog,Pdim,Pdeg,opts);
% Enforce strict positivity Pop>= eppos
eppos = 1e-3;
P = P +mat2opvar(eppos*eye(size(P)),P.dim,PIE.vars,PIE.dom);
% Indefinite operator variable Z
Zdim = B2.dim(:,[2,1]);
Zdeg = [4,0,0];
[prog,Z] = lpivar(prog,Zdim,Zdeg);

% % Set inequality constraints:
% %   Q <= 0
Iw = mat2opvar(eye(size(B1,2)), B1.dim(:,2), PIE.vars, PIE.dom);
Iz = mat2opvar(eye(size(C1,1)), C1.dim(:,1), PIE.vars, PIE.dom);
Q = [-gam*Iz,           D11,       (C1*P+D12*Z)*(T');
     D11',              -gam*Iw,   B1';
     T*(C1*P+D12*Z)',   B1,        (A*P+B2*Z)*(T')+T*(A*P+B2*Z)'];
prog = lpi_ineq(prog,-Q);

% % Set objective function:
% %   min gam
prog = lpisetobj(prog, gam);

% % Solve and retrieve the solution
opts.solver = 'sedumi';         % Use SeDuMi to solve the SDP
opts.simplify = true;           % Simplify the SDP before solving
prog_sol = lpisolve(prog,opts);
% Extract solved value of decision variables
gam_val = lpigetsol(prog_sol,gam);
Pval = lpigetsol(prog_sol,P);
Zval = lpigetsol(prog_sol,Z);
% Build the optimal control operator K.
Kval = getController(Pval,Zval,1e-3);

% % Construct the closed-loop PIE system
PIE_CL = pie_struct(); 
PIE_CL.vars = PIE.vars;         PIE_CL.dom = PIE.dom;
PIE_CL.T = T;
PIE_CL.A = A + B2*Kval;         PIE_CL.B1 = B1;
PIE_CL.C1 = C1 + D12*Kval;      PIE_CL.D11 = D11;
PIE_CL = initialize(PIE_CL);

% % Alternatively, uncomment to use pre-defined functions
% [prog, Kval, gam_val] = lpiscript(PIE,'hinf-controller','heavy');   
% PIE_CL = closedLoopPIE(PIE,Kval);


% =============================================
% === Simulate the system

% % Declare initial values and disturbance
syms st sx real
uinput1.ic.PDE = -10*sx;
uinput1.w = exp(-st);
uinput2.ic.PDE = sin(sx*pi/2);
uinput2.w = sin(pi*st)./(st+eps); 

% % Set options for discretization and simulation
opts.plot = 'no';   % don't plot final solution
opts.N = 16;        % expand using 16 Chebyshev polynomials
opts.tf = 5;        % simulate up to t = 5;
opts.dt = 1e-2;     % use time step of 10^-2
ndiff = [0,0,1];    % PDE state involves 1 second order differentiable state variable

% % Perform the actual simulation
% Simulate solution to uncontrolled PIE
[solution_OL_a,grid] = PIESIM(PIE,opts,uinput1,ndiff);
[solution_OL_b,~] = PIESIM(PIE,opts,uinput2,ndiff);
% Simulate solution to controlled PIE
[solution_CL_a,~] = PIESIM(PIE_CL,opts,uinput1,ndiff);
[solution_CL_b,~] = PIESIM(PIE_CL,opts,uinput2,ndiff);
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


echo off

% % Plot the open-loop state at different times
plot_indcs = [1,5,10,50,100];
colors = {'b','g','m','r','k','c','r','y','o'};     % colors for the plot
figure('Position',[200 150 1000 450]); 
for j=1:length(plot_indcs)
    plot_indx = plot_indcs(j);
    t_pos = num2str(tval(plot_indx));
    subplot(1,2,1)
    hold on
    plot(grid.phys,x_OL_a(:,plot_indx),[colors{j},'-'],'LineWidth',2,'DisplayName',['$\mathbf{x}_{ol}(t=',t_pos,')$']);
    hold off
    subplot(1,2,2)
    hold on
    plot(grid.phys,x_OL_b(:,plot_indx),[colors{j},'-'],'LineWidth',2,'DisplayName',['$\mathbf{x}_{ol}(t=',t_pos,')$']);
    hold off
end
% Clean up the figure
sgtitle('Open-loop PDE state','Interpreter','latex','FontSize',15)
subplot(1,2,1);
legend('Interpreter','latex','Location','northwest','FontSize',10.5);
xlabel('$s$','FontSize',15,'Interpreter','latex');  ylabel('$\mathbf{x}_{ol}$','FontSize',15,'Interpreter','latex');
title('$\mathbf{x}(0,s)=-\frac{5}{3}s^3+5s$,\quad $w(t)=e^{-t}$','Interpreter','latex','FontSize',13);
subplot(1,2,2);
legend('Interpreter','latex','Location','northwest','FontSize',10.5);
xlabel('$s$','FontSize',15,'Interpreter','latex');    ylabel('$\mathbf{x}_{ol}$','FontSize',15,'Interpreter','latex');
title('$\mathbf{x}(0,s)=-\frac{4}{\pi^2}\sin(\frac{\pi}{2}s)$,\quad $w(t)=\frac{\sin(\pi t)}{t}$','Interpreter','latex','FontSize',13);

% % Plot the closed-loop state at different times
figure('Position',[200 150 1000 450]); 
for j=1:length(plot_indcs)
    plot_indx = plot_indcs(j);
    t_pos = num2str(tval(plot_indx));
    subplot(1,2,1)
    hold on
    plot(grid.phys,x_CL_a(:,plot_indx),[colors{j},'-'],'LineWidth',2,'DisplayName',['$\mathbf{x}_{cl}(t=',t_pos,')$']);
    hold off
    subplot(1,2,2)
    hold on
    plot(grid.phys,x_CL_b(:,plot_indx),[colors{j},'-'],'LineWidth',2,'DisplayName',['$\mathbf{x}_{cl}(t=',t_pos,')$']);
    hold off
end
% Clean up the figure
sgtitle('Closed-loop PDE state','Interpreter','latex','FontSize',15)
subplot(1,2,1);
legend('Interpreter','latex','Location','southwest','FontSize',10.5);
xlabel('$s$','FontSize',15,'Interpreter','latex');  ylabel('$\mathbf{x}_{cl}$','FontSize',15,'Interpreter','latex');
title('$\mathbf{x}(0,s)=-\frac{5}{3}s^3+5s$,\quad $w(t)=e^{-t}$','Interpreter','latex','FontSize',13);
subplot(1,2,2);
legend('Interpreter','latex','Location','northwest','FontSize',10.5);
xlabel('$s$','FontSize',15,'Interpreter','latex');    ylabel('$\mathbf{x}_{cl}$','FontSize',15,'Interpreter','latex');
title('$\mathbf{x}(0,s)=-\frac{4}{\pi^2}\sin(\frac{\pi}{2}s)$,\quad $w(t)=\frac{\sin(\pi t)}{t}$','Interpreter','latex','FontSize',13);

% % Plot the closed-loop regulated output versus time
plot_indcs = floor(logspace(0,log(opts.tf/opts.dt)/log(10),50));
tplot = tval(plot_indcs);
figure('Position',[200 150 1000 450]); 
subplot(2,2,1)
hold on
plot(tplot,z_CL_a(1,plot_indcs),[colors{6},'-'],'LineWidth',1.5,'DisplayName','$z_{cl}(t)$');
plot(tplot,w_a(plot_indcs),[colors{7},'-.'],'LineWidth',1.0,'DisplayName','$w(t)$');
hold off
subplot(2,2,2)
hold on
plot(tplot,z_CL_b(1,plot_indcs),[colors{6},'-'],'LineWidth',1.5,'DisplayName','$z_{cl}(t)$');
plot(tplot,w_b(plot_indcs),[colors{7},'-.'],'LineWidth',1.0,'DisplayName','$w(t)$');
hold off
subplot(2,2,3)
plot(tplot,u_CL_a(1,plot_indcs),[colors{8},'-'],'LineWidth',1.5,'DisplayName','$u_{cl}(t)$');
subplot(2,2,4)
plot(tplot,u_CL_b(1,plot_indcs),[colors{8},'-'],'LineWidth',1.5,'DisplayName','$u_{cl}(t)$');
% Clean up the figure
sgtitle('Closed-loop regulated output $z_{cl}(t)=\int_{0}^{1}\mathbf{x}_{cl}(t,s)ds+w(t)$ and input effort $u_{cl}=\mathcal{K}\mathbf{x}_{cl}(t)$','Interpreter','latex','FontSize',15)
ax1 = subplot(2,2,1);   ax1.XScale = 'log';  ax1.XLim = [opts.dt,opts.tf];
legend('Interpreter','latex','Location','southwest','FontSize',10.5);
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