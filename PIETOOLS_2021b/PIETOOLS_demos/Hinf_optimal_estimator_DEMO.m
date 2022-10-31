% This document illustrates how an Hinfty optimal estimator can be designed
% for a PIE, and how the estimated state can be simulated.
%
% We consider the system defined by:
% PDE:      \dot{v}(t,s) = (d^2/ds^2) v(t,s) + 5 v(t,s) + w(t),    s in [0,1];
% Outputs:       z(t)    = int_{0}^{1} v(t,s) ds + w(t);
%                y(t)    = int_{0}^{1} v(t,s) ds;
% BCs:                0  = v(t,0) = v(t,1);
%
% Letting x:=(d^2/ds^2)v, We derive an equivalent PIE of the form:
% [T \dot{x}](t,s) = [A x](t,s) + [B1 w](t,s);
%             z(t) = [C1 x](t)  + [D11 w](t);
%             y(t) = [C2 x](t)  + [D21 w](t);
%
% We design an estimator of the form:
% [T \dot{xhat}](t,s) = [A xhat](t,s) + [L(yhat-y)](t,s);
%             zhat(t) = [C1 xhat](t)
%             yhat(t) = [C2 xhat](t)
%
% Using this estimator, the errors e=(xhat-x) and ztilde = (zhat-z) satisfy
% [T \dot{e}](t,s) = [(A+L*C2) e](t,s) - [(B1+L*D21) w](t,s);
%        ztilde(t) = [C1 e](t)         - [D11 w](t);
%
% We wish to compute an operator L that minimizes the L2 gain from
% disturbances w to error ztilde in the output. This is achieved solving
% the LPI
%
% min_{gam,P,Z} gam
% s.t.  P>=0
%       [-gam*I,           -D11',  -(P*B1+Z*D21)'*T           ]=: Q <=0
%       [-D11,             -gam*I, C1                         ]
%       [-T'*(P*B1+Z*D21), C1',    (P*A+Z*C2)'*T+T'*(P*A+Z*C2)]
%
% Then, using L = P^{-1}*Z, the L2 gain satisfied 
% ||ztilde||_{L2}/||w||_{L2} <= gam
%
% We manually declare this LPI here, but it can also be solved using the
% "PIETOOLS_Hinf_estimator" executive file.
% We simulate the PIE state x and estimated state xhat using PIESIM.
%
%%
clc; clear; close all;
echo on
%%%%%%%%%%%%%%%%%% Start Code Snippet %%%%%%%%%%%%%%%%%%
%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %%
% % % Declare the PDE, and convert it to a PIE.
% Declare the PDE using command line parser
pvar s t
PDE = sys();
x = state('pde');   w = state('in');
y = state('out');   z = state('out');
eqs = [diff(x,t) - diff(x,s,2) - 5*x - w;
       z - int(x,s,[0,1]) - w;
       y - int(x,s,[0,1]);
       subs(x,s,0);
       subs(x,s,1)];
PDE = addequation(PDE,eqs);
PDE = setObserve(PDE,y);
display_PDE(PDE);

% Compute the associated PIE, and extract the operators.
PIE = convert(PDE,'pie');       PIE = PIE.params;
T = PIE.T;
A = PIE.A;      C1 = PIE.C1;    C2 = PIE.C2;
B1 = PIE.B1;    D11 = PIE.D11;  D21 = PIE.D21;



%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %%
% % % Compute an optimal observer operator L for the PIE.

use_executive = false;  % <-- set to true to use predefined executive
if use_executive
    % % Use the predefined Hinf estimator executive function.
    settings = lpisettings('heavy');
    [prog, Lval, gam_val] = PIETOOLS_Hinf_estimator(PIE, settings);    
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
    Pdeg = {6,[2,3,5],[2,3,5]};
    opts.sep = 1;
    [prog,P] = poslpivar(prog,Pdim,Pdom,Pdeg,opts);
    %eppos = 1e-6;
    %P.R.R0 = P.R.R0 + eppos*eye(size(P));

    % Declare the indefinite PI operator decision variable Z
    Zdim = C2.dim(:,[2,1]);
    Zdom = PIE.dom;
    Zdeg = [4,0,0];
    [prog,Z] = lpivar(prog,Zdim,Zdom,Zdeg);

    % Declare the LPI constraint Q<=0.
    nw = size(B1,2);    nz = size(C1,1);
    Q = [-gam*eye(nw),     -D11',        -(P*B1+Z*D21)'*T;
         -D11,             -gam*eye(nz), C1;
         -T'*(P*B1+Z*D21), C1',          (P*A+Z*C2)'*T+T'*(P*A+Z*C2)];
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
    
    % Build the optimal observer operator L.
    Lval = getObserver(Pval,Zval);
end


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %%
% % % Build a PIE modeling the actual state x, and estimated state xhat:
% [T, 0] [\dot{x}(t)   ] = [A,     0     ] [x(t)   ] + [B1   ] w(t)
% [0, T] [\dot{xhat}(t)]   [-L*C2, A+L*C2] [xhat(t)] + [L*D21]
%
%              [z(t)   ] = [C1, 0 ] [x   ]           + [D11] w(t)
%              [zhat(t)]   [0,  C1] [xhat]             [0  ]

% Construct the operators defining the PIE.
T_CL = [T, 0*T; 0*T, T];
A_CL = [A, 0*A; -Lval*C2, A+Lval*C2];   B_CL = [B1; Lval*D21];
C_CL = [C1, 0*C1; 0*C1, C1];            D_CL = [D11; 0*D11];

% Declare the PIE.
PIE_CL = pie_struct();
PIE_CL.vars = PIE.vars;
PIE_CL.dom = PIE.dom;
PIE_CL.T = T_CL;        PIE_CL.Tw = opvar();    PIE_CL.Tu = opvar();
PIE_CL.A = A_CL;        PIE_CL.B1 = B_CL;       PIE_CL.B2 = opvar();
PIE_CL.C1 = C_CL;       PIE_CL.D11 = D_CL;      PIE_CL.D12 = opvar();
PIE_CL.C2 = opvar();    PIE_CL.D21 = opvar();   PIE_CL.D22 = opvar();


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %%
% % % Simulate and plot the actual and estimated output using PIESIM

% Declare initial conditions for the primary states of the PIE
syms st sx real
uinput.ic.PDE = [20*(sx^2-sx);  % Actual initial state value
                 0];            % Estimated initial state value
 
% Declare the value of the disturbance w(t)
uinput.w = 2*sin(pi*st);

% Set options for the discretization and simulation:
opts.plot = 'no';   % Do not plot the final solution
opts.N = 8;         % Expand using 8 Chebyshev polynomials
opts.tf = 1;        % Simulate up to t = 1;
opts.intScheme=1;   % Time-step using Backward Differentiation Formula (BDF) 
opts.dt = 1e-3;     % Use time step of 10^-3

% Simulate the solution to the PIE with estimator.
[solution,grid] = PIESIM(PIE_CL,opts,uinput,[0,0,2]);

% Extract actual and estimated solution at each time step.
x_act = reshape(solution.timedep.pde(:,1,:),opts.N+1,[]);
x_est = reshape(solution.timedep.pde(:,2,:),opts.N+1,[]);
tval = solution.timedep.dtime;

% Set options for the plot
plot_indcs = floor(linspace(1,opts.tf/opts.dt,76)); 
tplot = tval(plot_indcs);           % Only plot at select times
colors = {'g','b','r','m','k'};     % Colors for the plot
grid_idcs = [2,3,5];                % Only plot at a few grid points

% Plot evolution of actual and estimated
fig1 = figure(1);
hold on
for j=grid_idcs
    s_pos = num2str(grid.phys(j));  % Position associated to grid index.
    plot(tplot,x_act(j,plot_indcs),[colors{j},'-'],'LineWidth',2,'DisplayName',['$x(s=',s_pos,')$']);
    plot(tplot,x_est(j,plot_indcs),[colors{j},'o--'],'LineWidth',1.5,'DisplayName',['$\hat{x}(s=',s_pos,')$']);
end
hold off

% Clean up the figure
lgd = legend('Interpreter','latex');                lgd.FontSize=  10.5;
xlabel('Time','FontSize',15,'Interpreter','latex');  ylabel('State Value','FontSize',15,'Interpreter','latex');
title('Value of PIE state $x(t)$ and estimate $\hat{x}(t)$ at several grid points','Interpreter','latex','FontSize',16);
fig1.Position = [700 600 800 450];

%%
%%%%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
echo off