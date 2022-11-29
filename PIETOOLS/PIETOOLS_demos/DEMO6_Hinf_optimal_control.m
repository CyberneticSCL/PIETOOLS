% DEMO6_Hinf_optimal_control.m
% See Chapter 11.6 of the manual for a description.
%
% This document illustrates how an Hinfty optimal controller can be designed
% for a PDE using the PIE/LPI framework.
%
% We consider the system defined by:
% PDE:      \dot{x}(t,s) = (d^2/ds^2) x(t,s) + lam x(t,s) + w(t) + u(t),   s in [0,1];
% Outputs:       z(t)    = [int_{0}^{1} x(t,s) ds + w(t); u(t)];
% BCs:                0  = x(t,0) = d/ds x(t,1);
% (unstable for lam > pi^2/4)
%
% Letting v:=(d^2/ds^2)x, We derive an equivalent PIE of the form:
% [T \dot{v}](t,s) = [A v](t,s) + [B1 w](t,s) + [B2 u](t,s);
%             z(t) = [C1 v](t)  + [D11 w](t) + [D12 u](t);
%
% Using a state feedback control u = Kv we get the closed loop PIE
% [T \dot{v}](t,s) = [(A + B2*K) v](t,s) + [B1 w](t,s);
%             z(t) = [(C1 + D12*K) v](t) + [D11 w](t)
%
% We wish to compute an operator K that minimizes the L2 gain from
% disturbances w to the output z. This is achieved solving
% the LPI
%
% min_{gam,P,Z} gam
% s.t.  P>=0
%       [ -gam*I            D11      (C1*P+D12*Z)*T'            ]
%       [  D11'            -gam*I    B1'                        ] <= 0
%       [  T*(C1*P+D12*Z)   B1       (A*P+B2*Z)*T'+T*(A*P+B2*Z)']
%
% If the above LPI is successfully solved, then for K = Z*P^{-1}, we have
% ||z||_{L2}/||w||_{L2} <= gam
%
% We manually declare this LPI here, but it can also be solved using the
% "PIETOOLS_Hinf_control" executive file.
% We simulate the open loop and closed loop response of the PDE for various
% IC using PIESIM.
%
%%
clc; clear; close all;
echo on
%%%%%%%%%%%%%%%%%% Start Code Snippet %%%%%%%%%%%%%%%%%%
%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %%
% % % Declare the PDE, and convert it to a PIE.
% Declare the PDE using command line parser
pvar s t
lam = 5;
PDE = sys();
x = state('pde');   w = state('in');
z = state('out', 2);   u = state('in');
eqs = [diff(x,t) == diff(x,s,2) + lam*x + s*w + s*u;
    z == [int(x,s,[0,1]); u];
    subs(x,s,0)==0;
    subs(diff(x,s),s,1)==0];
PDE = addequation(PDE,eqs);
PDE = setControl(PDE,u);
display_PDE(PDE);

% Compute the associated PIE, and extract the operators.
PIE = convert(PDE,'pie');       PIE = PIE.params;
T = PIE.T;
A = PIE.A;      C1 = PIE.C1;    B2 = PIE.B2;
B1 = PIE.B1;    D11 = PIE.D11;  D12 = PIE.D12;



%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %%
% % % Compute an optimal observer operator L for the PIE.

use_executive = false;  % <-- set to true to use predefined executive
if use_executive
    % % Use the predefined Hinf estimator executive function.
    settings = lpisettings('heavy');
    [prog, Kval, gam_val] = PIETOOLS_Hinf_control(PIE, settings);
else
    % % Manually construct and solve the LPI program for optimal
    % % estimator synthesis.

    % Initialize the LPI program
    vars = PIE.vars(:);
    prog = sosprogram(vars);

    % Declare the decision variable gamma
    dpvar gam;
    prog = sosdecvar(prog, gam);

    % Declare a positive semidefinite PI operator decision variable P>=0
    Pdim = T.dim(:,1);
    Pdom = PIE.dom;
    Pdeg = {2,[1,1,2],[1,1,2]};
    opts.sep = 0;
    [prog,P] = poslpivar(prog,Pdim,Pdom,Pdeg,opts);
    eppos = 1e-3;
    P.R.R0 = P.R.R0 + eppos*eye(size(P.R.R0));

    % Declare the indefinite PI operator decision variable Z
    Zdim = B2.dim(:,[2,1]);
    Zdom = PIE.dom;
    Zdeg = [4,0,0];
    [prog,Z] = lpivar(prog,Zdim,Zdom,Zdeg);

    % Declare the LPI constraint Q<=0.
    nw = size(B1,2);    nz = size(C1,1);
    Q = [-gam*eye(nz)    D11          (C1*P+D12*Z)*(T');
        D11'           -gam*eye(nw)  B1';
        T*(C1*P+D12*Z)' B1           (A*P+B2*Z)*(T')+T*(A*P+B2*Z)'];
    prog = lpi_ineq(prog,-Q);

    % Set the objective function: minimize gam
    prog = sossetobj(prog, gam);

    % Solve the optimization program
    opts.solver = 'sedumi';
    opts.simplify = true;
    prog_sol = sossolve(prog,opts);

    % Extract the solved value of gam and the operators P and Z
    gam_val = sosgetsol(prog_sol,gam);
    Pval = getsol_lpivar(prog_sol,P);
    Zval = getsol_lpivar(prog_sol,Z);

    % Build the optimal observer operator K.
    Kval = getController(Pval,Zval,1e-3);
end


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %%

% Construct the operators defining the PIE.
use_CL_function = false;

if use_CL_function
    PIE_CL = closedLoopPIE(PIE,Kval);
else
    T_CL = T;
    A_CL = A+B2*Kval;   B_CL = B1;
    C_CL = C1+D12*Kval; D_CL = D11;

    % Declare the PIE.
    PIE_CL = pie_struct();
    PIE_CL.vars = PIE.vars;
    PIE_CL.dom = PIE.dom;
    PIE_CL.T = T_CL;
    PIE_CL.A = A_CL;        PIE_CL.B1 = B_CL;
    PIE_CL.C1 = C_CL;       PIE_CL.D11 = D_CL;
    PIE_CL = initialize(PIE_CL);
end

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %%
% % % Simulate and plot the actual and estimated PDE state using PIESIM

% Declare initial conditions for the state components of the PIE
syms st sx real

% Set options for the discretization and simulation:
opts.plot = 'no';   % Do not plot the final solution
opts.N = 8;         % Expand using 8 Chebyshev polynomials
opts.tf = 10;        % Simulate up to t = 1;
opts.dt = 1e-3;     % Use time step of 10^-3
opts.intScheme=1;   % Time-step using Backward Differentiation Formula (BDF)
ndiff = [0,0,1];    % The PDE state involves 1 second order differentiable state variables

% Simulate the solution to the PIE without controller for different IC.
uinput.ic.PDE = -10*sx;    % IC PIE
uinput.w = exp(-st); % disturbance
[solution_OL_a,grid] = PIESIM(PIE,opts,uinput,ndiff);
uinput.ic.PDE = sin(sx*pi/2);    % IC PIE
[solution_OL_b,grid] = PIESIM(PIE,opts,uinput,ndiff);


% Simulate the solution to the PIE with controller for different IC and disturbance.
uinput.ic.PDE = -10*sx;    % IC PIE
uinput.w = exp(-st);   % disturbance
[solution_CL_a,grid] = PIESIM(PIE_CL,opts,uinput,ndiff);
uinput.ic.PDE = sin(sx*pi/2);    % IC PIE
uinput.w = exp(-st);   % disturbance
[solution_CL_b,grid] = PIESIM(PIE_CL,opts,uinput,ndiff);

uinput.ic.PDE = -10*sx;    % IC PIE
uinput.w = sin(pi*st)./(st+eps);   % disturbance
[solution_CL_aa,grid] = PIESIM(PIE_CL,opts,uinput,ndiff);
uinput.ic.PDE = sin(sx*pi/2);    % IC PIE
uinput.w = sin(pi*st)./(st+eps);   % disturbance
[solution_CL_ab,grid] = PIESIM(PIE_CL,opts,uinput,ndiff);

%%
% Extract actual solution at each time step.
tval = solution_OL_a.timedep.dtime;
x_OL_a = reshape(solution_OL_a.timedep.pde(:,1,:),opts.N+1,[]);
x_OL_b = reshape(solution_OL_b.timedep.pde(:,1,:),opts.N+1,[]);
x_CL_a = reshape(solution_CL_a.timedep.pde(:,1,:),opts.N+1,[]);
x_CL_b = reshape(solution_CL_b.timedep.pde(:,1,:),opts.N+1,[]);
x_CL_aa = reshape(solution_CL_aa.timedep.pde(:,1,:),opts.N+1,[]);
x_CL_ab = reshape(solution_CL_ab.timedep.pde(:,1,:),opts.N+1,[]);
XX = linspace(0,1,20);
%%
% Set options for the plot
plot_indcs = floor(logspace(0,log(opts.tf/opts.dt)/log(10),5));
%plot_indcs = floor(linspace(1,opts.tf/opts.dt,66));
tplot = tval(plot_indcs);           % Only plot at select times
colors = {'b','g','m','r','k','c','r','y','o'};     % Colors for the plot
grid_idcs = 1:2:9;                % Only plot at a few grid points

% Plot open loop response for different IC
fig1 = figure(1); 
subplot(1,2,1); hold on;
for j = 1:length(plot_indcs)
    s_pos = num2str(plot_indcs(j)*opts.dt);  % Position associated to grid index.
    [YY] = spline(grid.phys,x_OL_a(:,plot_indcs(j)),XX);
    plot(XX,YY,[colors{j},'--o'],'LineWidth',2,'DisplayName',['$\mathbf{x}(t=',s_pos,')$'],'MarkerIndices',1:3:length(XX));
end 
hold off;
subplot(1,2,2); hold on;
for j = 1:length(plot_indcs)
    s_pos = num2str(plot_indcs(j)*opts.dt);  % Position associated to grid index.
    [YY] = spline(grid.phys,x_OL_b(:,plot_indcs(j)),XX);
    plot(XX,YY,[colors{j},'--o'],'LineWidth',2,'DisplayName',['$\mathbf{x}(t=',s_pos,')$'],'MarkerIndices',1:3:length(XX));
end
hold off;
fig1.Position = [100 100 3000 2000];
ax1 = subplot(1,2,1);     
set(ax1,'XTick',XX(1:4:end));
lgd1 = legend('Interpreter','latex');                lgd1.FontSize = 10.5;
lgd1.Location = 'northwest';
xlabel('$s$','FontSize',15,'Interpreter','latex');  ylabel('$\mathbf{x}(t,s)$','FontSize',15,'Interpreter','latex');
title('Open loop $\mathbf{x}(t)$; $\mathbf{x}(0)=\frac{5}{4}(1-s^2)$','Interpreter','latex','FontSize',15);
ax2 = subplot(1,2,2);
set(ax2,'XTick',XX(1:4:end));
lgd1 = legend('Interpreter','latex');                lgd2.FontSize = 10.5; 
lgd1.Location = 'northeast';
xlabel('$s$','FontSize',15,'Interpreter','latex');    ylabel('$\mathbf{x}(t,s)$','FontSize',15,'Interpreter','latex');
title('Open loop $\mathbf{x}(t)$; $\mathbf{x}(0)=-\frac{4}{\pi^2}\sin(\frac{\pi}{2} s)$','Interpreter','latex','FontSize',15);

%%
% Plot closed loop response for different IC
fig2 = figure(2); 
fig2.Position = [50 50 3000 2000];
subplot(2,2,1); hold on;
for j = 1:length(plot_indcs)
    s_pos = num2str(plot_indcs(j)*opts.dt);  % Position associated to grid index.
    [YY] = spline(grid.phys,x_CL_a(:,plot_indcs(j)),XX);
    plot(XX,YY,[colors{j},'--o'],'LineWidth',2,'DisplayName',['$\mathbf{x}(t=',s_pos,')$'],'MarkerIndices',1:3:length(XX));
end 
hold off;
subplot(2,2,2); hold on;
for j = 1:length(plot_indcs)
    s_pos = num2str(plot_indcs(j)*opts.dt);  % Position associated to grid index.
    [YY] = spline(grid.phys,x_CL_b(:,plot_indcs(j)),XX);
    plot(XX,YY,[colors{j},'--o'],'LineWidth',2,'DisplayName',['$\mathbf{x}(t=',s_pos,')$'],'MarkerIndices',1:3:length(XX));
end
hold off;
subplot(2,2,3); hold on;
for j = 1:length(plot_indcs)
    s_pos = num2str(plot_indcs(j)*opts.dt);  % Position associated to grid index.
    [YY] = spline(grid.phys,x_CL_aa(:,plot_indcs(j)),XX);
    plot(XX,YY,[colors{j},'--o'],'LineWidth',2,'DisplayName',['$\mathbf{x}(t=',s_pos,')$'],'MarkerIndices',1:3:length(XX));
end
hold off;
subplot(2,2,4); hold on;
for j = 1:length(plot_indcs)
    s_pos = num2str(plot_indcs(j)*opts.dt);  % Position associated to grid index.
    [YY] = spline(grid.phys,x_CL_ab(:,plot_indcs(j)),XX);
    plot(XX,YY,[colors{j},'--o'],'LineWidth',2,'DisplayName',['$\mathbf{x}(t=',s_pos,')$'],'MarkerIndices',1:3:length(XX));
end
hold off;

ax1 = subplot(2,2,1);
set(ax1,'XTick',XX(1:4:end));
lgd1 = legend('Interpreter','latex');                lgd1.FontSize = 10.5;
lgd1.Location = 'southwest';
xlabel('$s$','FontSize',15,'Interpreter','latex');  ylabel('$\mathbf{x}(t,s)$','FontSize',15,'Interpreter','latex');
title('Closed loop $\mathbf{x}(t)$; $\mathbf{x}(0)=\frac{5}{4}(1-s^2)$, $w(t)=e^{-t}$','Interpreter','latex','FontSize',15);
ax2 = subplot(2,2,2);
set(ax2,'XTick',XX(1:4:end));
lgd1 = legend('Interpreter','latex');                lgd2.FontSize = 10.5; 
lgd1.Location = 'northwest';
xlabel('$s$','FontSize',15,'Interpreter','latex');    ylabel('$\mathbf{x}(t,s)$','FontSize',15,'Interpreter','latex');
title('Closed loop $\mathbf{x}(t)$; $\mathbf{x}(0)=-\frac{4}{\pi^2}\sin(\frac{\pi}{2}s)$, $\mathbf{x}(0)=e^{-t}$','Interpreter','latex','FontSize',15);
ax3 = subplot(2,2,3);
set(ax3,'XTick',XX(1:4:end));
lgd1 = legend('Interpreter','latex');                lgd1.FontSize = 10.5;
lgd1.Location = 'southwest';
xlabel('$s$','FontSize',15,'Interpreter','latex');  ylabel('$\mathbf{x}(t,s)$','FontSize',15,'Interpreter','latex');
title('Closed loop $\mathbf{x}(t)$; $\mathbf{x}(0)=\frac{5}{4}(1-s^2)$, $w(t)=\frac{\sin(\pi t)}{t}$','Interpreter','latex','FontSize',15);
ax4 = subplot(2,2,4);     
set(ax4,'XTick',XX(1:4:end));
lgd1 = legend('Interpreter','latex');                lgd2.FontSize = 10.5; 
lgd1.Location = 'northwest';
xlabel('$s$','FontSize',15,'Interpreter','latex');    ylabel('$\mathbf{x}(t,s)$','FontSize',15,'Interpreter','latex');
title('Closed loop $\mathbf{x}(t)$; $\mathbf{x}(0)=-\frac{4}{\pi^2}\sin(\frac{\pi}{2}s)$, $w(t)=\frac{\sin(\pi t)}{t}$','Interpreter','latex','FontSize',15);
%%
w1_tval = subs(sin(pi*st)./(st+eps),tval);
w2_tval = subs(exp(-st),tval);
z_quadrature = double(subs(0.5*sx^2-sx,grid.phys));
k_quadrature = double(subs(Kval.Q1,s,grid.phys));
ZZ1 = trapz(z_quadrature,x_CL_a)+double(w1_tval); %int wa
ZZ2 = trapz(z_quadrature,x_CL_aa)+double(w2_tval);%int wb
ZZ3 = trapz(k_quadrature,x_CL_a); %u wa
ZZ4 = trapz(k_quadrature,x_CL_aa); %u wb
fig3 = figure(3); XX = linspace(0,1,2000);
subplot(1,2,1); hold on;
[YY] = spline(tval,ZZ1,XX);
plot(XX,YY,'--o','MarkerIndices',1:90:length(XX),'LineWidth',2,'DisplayName',['$w(t) = e^{-t}$']);
[YY] = spline(tval,ZZ2,XX);
plot(XX,YY,'--x','MarkerIndices',1:90:length(XX),'LineWidth',2,'DisplayName',['$w(t) = \frac{\sin(\pi t)}{t}$']); hold off;
subplot(1,2,2); hold on;
[YY] = spline(tval,ZZ3,XX);
plot(XX,YY,'--o','MarkerIndices',1:90:length(XX),'LineWidth',2,'DisplayName',['$w(t) = e^{-t}$']);
[YY] = spline(tval,ZZ4,XX);
plot(XX,YY,'--x','MarkerIndices',1:90:length(XX),'LineWidth',2,'DisplayName',['$w(t) = \frac{\sin(\pi t)}{t}$']); hold off;

ax5 = subplot(1,2,1);
set(ax5,'XTick',tval(1:90:end),'xlim',[0,0.75]);
lgd1 = legend('Interpreter','latex'); lgd1.FontSize = 10.5; 
lgd1.Location = 'southeast';
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$z(t)$','FontSize',15,'Interpreter','latex');
title('Closed loop $z(t) = \int_0^1 \mathbf{x}(t,s) ds+w(t)$','Interpreter','latex','FontSize',15);
ax6 = subplot(1,2,2);
set(ax6,'XTick',tval(1:90:end),'xlim',[0,0.75]);
lgd1 = legend('Interpreter','latex'); lgd1.FontSize = 10.5; 
lgd1.Location = 'southeast';
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$u(t)$','FontSize',15,'Interpreter','latex');
title('Control effort $u(t)$','Interpreter','latex','FontSize',15);
%%
%%%%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
echo off
