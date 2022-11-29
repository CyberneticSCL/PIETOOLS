% DEMO7_observer_based_control.m
% See Chapter 11.7 of the manual for a description.
%
% This document illustrates construction of an observer-based state feedback 
% controller. We construct the observer and controller separately
% and couple them to find a closed loop system which is then simulated for 
% different IC and disturbances.
%
% We consider the wave equation PDE defined by:
% PDE:      \dot{x}(t,s) = (d^2/ds^2) x(t,s) + lam x(t,s) + w(t) + u(t),   s in [0,1];
% Outputs:       z(t)    = [int_{0}^{1} x(t,s) ds + w(t); u(t)];
%                y(t0    = x(t,1);
% BCs:                0  = x(t,0) = d/ds x(t,1);
%
% First we convert the above PDE to a PIE of the form:
% [T \dot{v}](t,s) = [A v](t,s) + [B1 w](t,s) + [B2 u](t,s);
%             z(t) = [C1 v](t)  + [D11 w](t) + [D12 u](t);
%             y(t) = [C2 v](t)  + [D21 w](t) + [D22 u](t)
%
% We find the observer gains L by solving the LPI
%
% min_{gam,P,Z} gam
% s.t.  P>=0
%       [-gam*I,           -D11',  -(P*B1+Z*D21)'*T           ]=: Q <=0
%       [-D11,             -gam*I, C1                         ]
%       [-T'*(P*B1+Z*D21), C1',    (P*A+Z*C2)'*T+T'*(P*A+Z*C2)]
%
% where L = P^{-1}*Z.

% Likewise we find a state feedback control u = Kv by solving
% the LPI
%
% min_{gam,P,Z} gam
% s.t.  P>=0
%       [ -gam*I            D11      (C1*P+D12*Z)*T'            ]
%       [  D11'            -gam*I    B1'                        ] <= 0
%       [  T*(C1*P+D12*Z)   B1       (A*P+B2*Z)*T'+T*(A*P+B2*Z)']
%
% where K = Z*P^{-1}.
%
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
x = state('pde');   w = state('in'); y = state('out');
z = state('out', 2);   u = state('in');
eqs = [diff(x,t) == diff(x,s,2) + lam*x + s*w + s*u;
    z == [int(x,s,[0,1]); u]; % change z = int(x,s,[0,1])
    y == subs(x,s,1);
    subs(x,s,0)==0;
    subs(diff(x,s),s,1)==0];
PDE = addequation(PDE,eqs);
PDE = setControl(PDE,u);
PDE = setObserve(PDE,y);
display_PDE(PDE);

% Compute the associated PIE, and extract the operators.
PIE = convert(PDE,'pie');       PIE = PIE.params;
T = PIE.T;
A = PIE.A;      C1 = PIE.C1;    B2 = PIE.B2;
B1 = PIE.B1;    D11 = PIE.D11;  D12 = PIE.D12;
C2 = PIE.C2;    D21 = PIE.D21;  D22 = PIE.D22;


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %%
% % % Compute an optimal observer operator L for the PIE.

use_executive = true;  % <-- set to true to use predefined executive
if use_executive
    % % Use the predefined Hinf estimator executive function.
    settings = lpisettings('heavy');
    [prog_k, Kval, gam_co_val] = PIETOOLS_Hinf_control(PIE, settings);
    [prog_l, Lval, gam_ob_val] = PIETOOLS_Hinf_estimator(PIE, settings);
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
    P.R.R0 = P.R.R0 + eppos*eye(size(P));

    % Declare the indefinite PI operator decision variable Z
    Zdim = B2.dim(:,[2,1]);
    Zdom = PIE.dom;
    Zdeg = [4,0,0];
    [prog,Z_con] = lpivar(prog,Zdim,Zdom,Zdeg);
    
    % Declare the LPI constraint Q_con<=0 for controller.
    nw = size(B1,2);    nz = size(C1,1);
    Q_con = [-gam*eye(nz)    D11          (C1*P+D12*Z_con)*(T');
        D11'           -gam*eye(nw)  B1';
        T*(C1*P+D12*Z_con)' B1           (A*P+B2*Z_con)*(T')+T*(A*P+B2*Z_con)'];
    prog = lpi_ineq(prog,-Q_con);

    Zdim = C2.dim(:,[2,1]);
    Zdom = PIE.dom;
    Zdeg = [4,0,0];
    [prog,Z_obs] = lpivar(prog,Zdim,Zdom,Zdeg);

    % Declare the LPI constraint Q_obs<=0 for observer.
    nw = size(B1,2);    nz = size(C1,1);
    Q_obs = [-gam*eye(nw),     -D11',        -(P*B1+Z_obs*D21)'*T;
         -D11,             -gam*eye(nz), C1;
         -T'*(P*B1+Z_obs*D21), C1',          (P*A+Z_obs*C2)'*T+T'*(P*A+Z_obs*C2)];
    prog = lpi_ineq(prog,-Q_obs);


    % Set the objective function: minimize gam
    prog = sossetobj(prog, gam);

    % Solve the optimization program
    opts.solver = 'sedumi';
    opts.simplify = true;
    prog_sol = sossolve(prog,opts);

    % Extract the solved value of gam and the operators P and Z
    gam_val = sosgetsol(prog_sol,gam);
    Pval = getsol_lpivar(prog_sol,P);
    Zval_con = getsol_lpivar(prog_sol,Z_con);
    Zval_obs = getsol_lpivar(prog_sol,Z_obs);

    % Build the optimal observer operator K.
    Kval = getController(Pval,Zval_con,1e-3);
    Lval = getObserver(Pval,Zval_obs,1e-3);
end


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %%

% Construct the operators defining the PIE.
T_CL = [T, 0*T; 0*T, T];
A_CL = [A, B2*Kval; -Lval*C2, A+Lval*C2];   B_CL = [B1; Lval*D21];
C_CL = [C1, D12*Kval; 0*C1, C1];            D_CL = [D11; 0*D11];



% Declare the PIE.
PIE_CL = pie_struct();
PIE_CL.vars = PIE.vars;
PIE_CL.dom = PIE.dom;
PIE_CL.T = T_CL;
PIE_CL.A = A_CL;        PIE_CL.B1 = B_CL;
PIE_CL.C1 = C_CL;       PIE_CL.D11 = D_CL;
PIE_CL = initialize(PIE_CL);

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
uinput.ic.PDE = [-10*sx];   % IC PIE 
uinput.w = 5*exp(-st); % disturbance
[solution_OL,grid] = PIESIM(PIE,opts,uinput,ndiff);

% Simulate the solution to the PIE with controller for different IC and disturbance.
ndiff = [0,0,2]; 
uinput.ic.PDE = [-10*sx; 0];    % IC PIE and observed state
uinput.w = 5*exp(-st);   % disturbance
[solution_CL_a,grid] = PIESIM(PIE_CL,opts,uinput,ndiff);
uinput.ic.PDE = [sin(sx*pi/2); 0];    % IC PIE
uinput.w = 5*sin(pi*st)./(st+eps);   % disturbance
[solution_CL_b,grid] = PIESIM(PIE_CL,opts,uinput,ndiff);


%%
% Extract actual solution at each time step.
tval = solution_OL.timedep.dtime;
x_OL = reshape(solution_OL.timedep.pde(:,1,:),opts.N+1,[]);
x_CL_a = reshape(solution_CL_a.timedep.pde(:,1,:),opts.N+1,[]);
hatx_CL_a = reshape(solution_CL_a.timedep.pde(:,2,:),opts.N+1,[]);
x_CL_b = reshape(solution_CL_b.timedep.pde(:,1,:),opts.N+1,[]);
hatx_CL_b = reshape(solution_CL_b.timedep.pde(:,2,:),opts.N+1,[]);
%%
% Set options for the plot
plot_indcs = floor(logspace(0,log(opts.tf/opts.dt)/log(10),5));
%plot_indcs = floor(linspace(1,opts.tf/opts.dt,66));
tplot = tval(plot_indcs);           % Only plot at select times
colors = {'b','g','m','r','k','c','r','y','o'};     % Colors for the plot
grid_idcs = 1:2:9;                % Only plot at a few grid points

%%
tsteps = 1:50:length(tval);
fig1 = figure(1);
surf(tval(tsteps),grid.phys,x_OL(:,tsteps,:),'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
h=colorbar ;
colormap jet
box on
ylabel(h,'$|\mathbf{x}(t,s)|$','interpreter', 'latex','FontSize',15)
set(gcf, 'Color', 'w');
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$s$','FontSize',15,'Interpreter','latex');
zlabel('$\mathbf{x}(t,s)$','FontSize',15,'Interpreter','latex');
title('Open loop response with $\mathbf{x}(0,s)=\frac{5}{4}(1-s^2)$, $w(t)=5e^{-t}$','Interpreter','latex','FontSize',15);

fig2 = figure(2);
subplot(2,1,1);
surf(tval(tsteps),grid.phys,x_CL_a(:,tsteps,:),'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
h=colorbar ;
colormap jet
box on
ylabel(h,'$|\mathbf{x}(t,s)|$','interpreter', 'latex','FontSize',15)
set(gcf, 'Color', 'w'); 
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$s$','FontSize',15,'Interpreter','latex');
zlabel('$\mathbf{x}(t,s)$','FontSize',15,'Interpreter','latex');
title('Closed loop response with $\mathbf{x}(0,s)=\frac{5}{4}(1-s^2)$, $w(t)=5e^{-t}$','Interpreter','latex','FontSize',15);

subplot(2,1,2);
surf(tval(tsteps),grid.phys,x_CL_b(:,tsteps,:),'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
h=colorbar ;
colormap jet
box on
ylabel(h,'$|\mathbf{x}(t,s)|$','interpreter', 'latex','FontSize',15)
set(gcf, 'Color', 'w');
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$s$','FontSize',15,'Interpreter','latex');
zlabel('$\mathbf{x}(t,s)$','FontSize',15,'Interpreter','latex');
title('Closer loop response with $\mathbf{x}(0,s)=-\frac{4}{\pi^2}\sin(\frac{\pi}{2}s)$, $w(t)=5\frac{\sin(\pi t)}{t}$','Interpreter','latex','FontSize',15);
%%
XX = linspace(0,10,200);
w2_tval = subs(5*sin(pi*st)./(st+eps),tval);
w1_tval = subs(5*exp(-st),tval);
z_quadrature = double(subs(0.5*sx^2-sx,grid.phys));
k_quadrature = double(subs(Kval.Q1,s,grid.phys));
ZZ1 = trapz(z_quadrature,x_CL_a); %int wa
ZZ2 = trapz(z_quadrature,x_CL_b);%int wb
ZZ3 = trapz(k_quadrature,hatx_CL_a); %u wa
ZZ4 = trapz(k_quadrature,hatx_CL_b); %u wb
fig3 = figure(3); 
subplot(1,2,1); hold on;
[YY] = spline(tval,ZZ1,XX);
plot(XX,YY,'--o','MarkerIndices',1:90:length(XX),'LineWidth',2,'DisplayName',['$w(t) = 5e^{-t}$']);
[YY] = spline(tval,ZZ2,XX);
plot(XX,YY,'--x','MarkerIndices',1:90:length(XX),'LineWidth',2,'DisplayName',['$w(t) = 5\frac{\sin(\pi t)}{t}$']); hold off;
subplot(1,2,2); hold on;
[YY] = spline(tval,ZZ3,XX);
plot(XX,YY,'--o','MarkerIndices',1:90:length(XX),'LineWidth',2,'DisplayName',['$w(t) = e^{-t}$']);
[YY] = spline(tval,ZZ4,XX);
plot(XX,YY,'--x','MarkerIndices',1:90:length(XX),'LineWidth',2,'DisplayName',['$w(t) = 5\frac{\sin(\pi t)}{t}$']); hold off;

ax5 = subplot(1,2,1);
set(ax5,'XTick',tval(1:900:end));
lgd1 = legend('Interpreter','latex'); lgd1.FontSize = 10.5; 
lgd1.Location = 'northeast';
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$z(t)$','FontSize',15,'Interpreter','latex');
title('Closed loop $z(t) = \int_0^1 \mathbf{x}(t,s) ds$','Interpreter','latex','FontSize',15);
ax6 = subplot(1,2,2);
set(ax6,'XTick',tval(1:900:end));
lgd1 = legend('Interpreter','latex'); lgd1.FontSize = 10.5; 
lgd1.Location = 'southeast';
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$u(t)$','FontSize',15,'Interpreter','latex');
title('Control effort $u(t)$','Interpreter','latex','FontSize',15);
%%
%%%%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
echo off
