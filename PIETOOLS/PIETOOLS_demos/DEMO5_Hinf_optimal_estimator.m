% DEMO5_Hinf_optimal_estimator.m
% See Chapter 11.5 of the manual for a description.
%
% This document illustrates how an Hinfty optimal estimator can be designed
% for a PIE, and how the estimated state can be simulated.
% Specifically, we consider the following system:
% PDE:      \dot{x}(t,s) = (d^2/ds^2) x(t,s) + lam x(t,s) + w(t),   s in [0,1];
% Outputs:          z(t) = int_{0}^{1} x(t,s) ds + w(t);
%                   y(t) = x(t,1);
% BCs:                0  = x(t,0) = x_{s}(t,1);
% (unstable for lam > pi^2/4 = 2.4674)
% Letting v:=(d^2/ds^2)x, We derive an equivalent PIE of the form:
%   [T \dot{v}](t,s) = [A v](t,s) + [B1 w](t,s);
%               z(t) = [C1 v](t)  + [D11 w](t);
%               y(t) = [C2 v](t)  + [D21 w](t);
% So that x = T*v. We design an estimator of the form:
%   [T \dot{vhat}](t,s) = [A vhat](t,s) + [L(yhat-y)](t,s);
%               zhat(t) = [C1 vhat](t)
%               yhat(t) = [C2 vhat](t)
% Then, the errors e=(vhat-v) and ztilde = (zhat-z) satisfy
%   [T \dot{e}](t,s) = [(A+L*C2) e](t,s) - [(B1+L*D21) w](t,s);
%          ztilde(t) = [C1 e](t)         - [D11 w](t);
% To compute an operator L that minimizes the L2 gain from
% disturbances w to error ztilde in the output, we solve the LPI
%   min_{gam,P,Z}   gam
%   s.t.            P>=0,
%       [-gam*I,           -D11',  -(P*B1+Z*D21)'*T           ]=: Q <=0
%       [-D11,             -gam*I, C1                         ]
%       [-T'*(P*B1+Z*D21), C1',    (P*A+Z*C2)'*T+T'*(P*A+Z*C2)]
%
% Then, using L = P^{-1}*Z, the L2 gain satisfies 
%   ||ztilde||_{L2}/||w||_{L2} <= gam
% We manually declare this LPI here, but it can also be solved using the
% "PIETOOLS_Hinf_estimator" executive file.
% We simulate the PDE state x and estimated state xhat using PIESIM.
%

clc; clear; close all;
echo on

% =============================================
% === Declare the system of interest

% % Declare system as PDE
% Declare independent variables (time and space)
pvar s t
% Declare state, input, and output variables
x = state('pde');   w = state('in');
y = state('out');   z = state('out');
% Declare the sytem equations
lam = 4;
PDE = sys();
eqs = [diff(x,t) == diff(x,s,2) + lam*x + w;    % PDE
       z == int(x,s,[0,1]) + w;                 % regulated output
       y == subs(x,s,1);                        % observed output
       subs(x,s,0) == 0;                        % first boundary condition
       subs(diff(x,s),s,1) == 0];               % second boundary condition
PDE = addequation(PDE,eqs);                     % Add the equations to the PDE structure
% Set y as an observed output
PDE = setObserve(PDE,y);                        
display_PDE(PDE);

% % Convert PDE to PIE
PIE = convert(PDE);
T = PIE.T;
A = PIE.A;      C1 = PIE.C1;    C2 = PIE.C2;
B1 = PIE.B1;    D11 = PIE.D11;  D21 = PIE.D21;


% =============================================
% === Declare the LPI

% % Initialize LPI program
prog = lpiprogram(PIE.vars,PIE.dom);

% % Declare decision variables:
% %   gam \in \R,     P:L2-->L2,    Z:\R-->L2
% Scalar decision variable
[prog,gam] = lpidecvar(prog,'gam');
% Positive operator variable P>=0
opts.sep = 1;                   % set P.R.R1=P.R.R2
[prog,P] = poslpivar(prog,T.dim,4,opts);
% Indefinite operator variable Z: \R-->L2
[prog,Z] = lpivar(prog,[0,1;1,0],4);

% % Set inequality constraints:
% %   Q <= 0
Q = [-gam,              -D11',        -(P*B1+Z*D21)'*T;
     -D11,              -gam,         C1;
     -T'*(P*B1+Z*D21),  C1',          (P*A+Z*C2)'*T+T'*(P*A+Z*C2)];
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
% Build optimal observer operator L
Lval = getObserver(Pval,Zval);

% % Construct the closed-loop PIE system
% % T_CL * \dot{V}(t) = A_CL * V(t) + B_CL * w(t)
%               Z(t)  = C_CL * V(t) + D_CL * w(t)
% % where V = [v; vhat] and Z = [z; zhat].
PIE_CL = pie_struct(); 
PIE_CL.vars = PIE.vars;         PIE_CL.dom = PIE.dom;                       
PIE_CL.T = [T, 0*T; 0*T, T];
PIE_CL.A = [A, 0*A; -Lval*C2, A+Lval*C2];
PIE_CL.B1 = [B1; Lval*D21];
PIE_CL.C1 = [C1, 0*C1; 0*C1, C1];       
PIE_CL.D11 = [D11; 0*D11]; 
PIE_CL = initialize(PIE_CL);

% % Alternatively, uncomment to use pre-defined functions
% [prog, Lval, gam_val] = lpiscript(PIE,'hinf-observer','heavy');   
% PIE_CL = closedLoopPIE(PIE,Lval,'observer');


% =============================================
% === Simulate the system

% % Declare initial values and disturbance
syms st sx real
uinput.ic.PDE = [-10*sx;    % actual initial PIE state value
                 0];        % estimated initial PIE state value
uinput.w = 2*sin(pi*st);

% % Set options for discretization and simulation
opts.plot = 'no';   % don't plot final solution
opts.N = 8;         % expand using 8 Chebyshev polynomials
opts.tf = 1;        % simulate up to t = 1;
opts.dt = 1e-3;     % use time step of 10^-3
opts.intScheme = 1; % time-step using Backward Differentiation Formula (BDF) 
ndiff = [0,0,2];    % PDE state involves 2 second order differentiable state variables

% % Simulate and plot solution to the PIE with estimator.
[solution,grid] = PIESIM(PIE_CL,opts,uinput,ndiff);
% Extract actual and estimated state at each time step.
tval = solution.timedep.dtime;
x_act = reshape(solution.timedep.pde(:,1,:),opts.N+1,[]);
x_est = reshape(solution.timedep.pde(:,2,:),opts.N+1,[]);


echo off


% % Plot actual and estimated state at select positions and times
grid_idcs = [1,5,7];        % grid points at which to plot
plot_indcs = floor(logspace(0,log(opts.tf/opts.dt)/log(10),40)); 
tplot = tval(plot_indcs);   % times at which to plot
colors = {'b','g','m','r','k','c','r','y','o'};
figure('Position',[200 150 1000 450]);
for j=grid_idcs
    s_pos = num2str(grid.phys(j));  % position associated to grid index
    subplot(1,2,1)
    hold on
    plot(tplot,x_act(j,plot_indcs),[colors{j},'-'],'LineWidth',2,'DisplayName',['$\mathbf{x}(s=',s_pos,')$']);
    plot(tplot,x_est(j,plot_indcs),[colors{j},'o--'],'LineWidth',1.5,'DisplayName',['$\mathbf{\hat{x}}(s=',s_pos,')$']);
    hold off
    subplot(1,2,2)
    hold on
    loglog(tplot,(x_act(j,plot_indcs)-x_est(j,plot_indcs)),[colors{j},'o-'],'LineWidth',1.5,'DisplayName',['$\mathbf{e}(s=',s_pos,')$']);
    hold off
end
% Clean up the figure
ax1 = subplot(1,2,1);   ax1.XScale = 'log';     ax1.XTickLabels = {'0.001';'0.01';'0.1';'1'};
legend('Interpreter','latex','Location','northwest','FontSize',10.5);
xlabel('$t$','FontSize',15,'Interpreter','latex');  ylabel('$\mathbf{x}$','FontSize',15,'Interpreter','latex');
title('PDE state value $\mathbf{x}(t)$ and estimate $\mathbf{\hat{x}}(t)$','Interpreter','latex','FontSize',15);
ax2 = subplot(1,2,2);   ax2.XScale = 'log';     ax2.XTickLabels = {'0.001';'0.01';'0.1';'1'};
legend('Interpreter','latex','FontSize',10.5);
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$\mathbf{e}$','FontSize',15,'Interpreter','latex');
title('Error $\mathbf{e}=\mathbf{\hat{x}}(t)-\mathbf{x}(t)$','Interpreter','latex','FontSize',15);