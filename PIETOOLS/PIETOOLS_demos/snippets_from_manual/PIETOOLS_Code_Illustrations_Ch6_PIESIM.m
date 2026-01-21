% This document illustrates how PIESIM can be used to simulate solutions
% to 1D PDEs, 2D PDEs, DDEs and PIEs.
% We refer to Chapter 6 of the manual for more context on the codes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Code Illustrations
%
% Copyright (C)2024  PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, make sure to change the code in the manual as
% well, and vice versa. Document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 01/01/2025: Initial coding;
clc; clear; close all; clear stateNameGenerator

%% 6.4.1: 1D PDE Example
% Declare the PDE to simulate
[PDE,uinput] = examples_pde_library_PIESIM_1D(4);
% Declare initial conditions and disturbance
syms sx st;
uinput.exact(1) = -2*sx*st-sx^2;
uinput.w(1) = 0;        % disturbance at lower boundary
uinput.w(2) = -4*st-4;  % disturbance at upper boundary
uinput.ic.PDE = -sx^2;
% Declare discretization and temporal integration options
opts.ifexact = true;
opts.plot = 'yes';
opts.N = 8;
opts.tf = 1;
opts.Norder = 2;
opts.dt = 1e-3;
[solution,grids] = PIESIM(PDE,opts,uinput);
tval = solution.timedep.dtime;
xval = reshape(solution.timedep.pde(:,1,:),opts.N+1,[]);
% % Extract the solution
plot_indcs = floor(linspace(1,opts.tf/opts.dt,100));
tplot = tval(plot_indcs);
x_sim = xval(:,plot_indcs);
x_true = double(subs(subs(uinput.exact(1),st,tplot),sx,grids.phys(:)));

% Plot the evolution of the PDE state, as well as the final value.
fig1 = figure('Position',[200 50 900 350]); 
set(gcf, 'Color', 'w');
subplot(1,2,1);
box on
surf(grids.phys,tplot,x_sim','FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','column');
h = colorbar;
%colormap jet    
subplot(1,2,2);
box on
plot(grids.phys(:),x_sim(:,end),'rd','LineWidth',4,'DisplayName','Numerical solution');
hold on
plot(grids.phys(:),x_true(:,end),'k-','LineWidth',2,'DisplayName','Analytic solution');
hold off
%colormap jet
% Clean up the figure
sgtitle('Simulated PDE State Evolution','Interpreter','latex','FontSize',17)
ax1 = subplot(1,2,1);       ax1.TickLabelInterpreter = 'latex';
ax1.Position = [0.08,0.12,0.33,0.68];       ax1.ZLim = [-8,0];
box on
title('Evolution of PDE state $\mathbf{x}(t,s)$','Interpreter','latex','FontSize',15);
xlabel('$s$','FontSize',15,'Interpreter','latex');  
ylabel('$t$','FontSize',15,'Interpreter','latex');
zlabel('$\mathbf{x}$','FontSize',15,'Interpreter','latex');
ax2 = subplot(1,2,2);       ax2.TickLabelInterpreter = 'latex';
ax2.Position = [0.6,0.15,0.33,0.65];
lgd = legend('Location','northeast','Interpreter','latex','FontSize',15);
box on
title('Final PDE state $\mathbf{x}(t=1,s)$','Interpreter','latex','FontSize',15);
xlabel('$s$','FontSize',15,'Interpreter','latex');  
ylabel('$\mathbf{x}$','FontSize',15,'Interpreter','latex');
%saveas(fig1,'Ch6_ExA_1D_Plot','epsc');



%% 6.4.2: 2D PDE Example
clear PDE uinput
% Declare the desired 2D PDE
[PDE,~] = examples_pde_library_PIESIM_2D(9);

% % Set the exact solution
% %   u(x,y,t) = 10*cos(pi*(x-a)/(b-a))*cos(5*pi*(y-c)/(d-c))*exp(-t)
syms sx sy st real
u_ex = 10*cos(sym(pi)*sx)*cos(sym(pi)*sy)*exp(-st);
uinput.exact(1) = subs(subs(u_ex,sy,0),sx,0);
uinput.exact(2) = subs(u_ex,sy,0);
uinput.exact(3) = subs(u_ex,sx,0);
uinput.exact(4) = u_ex;
% % Set the initial conditions
uinput.ic.x = subs(uinput.exact,st,0);

% % Set the simulation options and simulate
opts.ifexact = true;
opts.plot = 'yes';
opts.N = 16;
opts.tf = 1;
opts.intScheme = 1;
opts.dt = 1e-3;
[solution,grids] = PIESIM(PDE,opts,uinput);
% Extract numerical solution at final time
x1val = solution.timedep.ode;         % x1(t)
x2fin = solution.final.pde{1}(:,1);   % x2(t=tf,s1);
x3fin = solution.final.pde{1}(:,2);   % x3(t=tf,s2);
x4fin = solution.final.pde{2};        % x4(t=tf,s1,s2);
% Extract grid points and exact solution
s1_grid = grids.phys(:,1);
s2_grid = grids.phys(:,2);
s1_grid_exact = linspace(0,1,round(100));
s2_grid_exact = linspace(0,5,round(100));
x2fin_true = double(subs(subs(uinput.exact(2),st,opts.tf),sx,s1_grid_exact));
x3fin_true = double(subs(subs(uinput.exact(3),st,opts.tf),sy,s2_grid_exact));
x4fin_true = double(subs(subs(subs(uinput.exact(4),st,opts.tf),sy,s2_grid_exact),sx,s1_grid_exact'));

% Plot 1D state values at final time
fig1 = figure('Position',[200 50 900 350]); 
sgtitle('1D PDE States at Final Time','FontSize',15,'Interpreter','latex');
set(gcf, 'Color', 'w');
ax1 = subplot(1,2,1);
box on
plot(s1_grid,x2fin,'rd','LineWidth',4,'DisplayName','Numerical solution');
hold on
plot(s1_grid_exact,x2fin_true,'k-','LineWidth',2,'DisplayName','Analytic solution');
hold off
legend('FontSize',13,'Interpreter','latex');
xlabel('$s_{1}$','FontSize',15,'Interpreter','latex');
ylabel('$\mathbf{x}_{2}$','FontSize',15,'Interpreter','latex');
title(['$\mathbf{x}_{2}(t=',num2str(solution.timedep.dtime(end)),',s_{1})$'],'FontSize',14,'Interpreter','latex');
set(gca,'XLim',[min(s1_grid),max(s1_grid)]);
set(gca,'TickLabelInterpreter','latex');
ax1.Position = [0.07,0.15,0.4,0.72];
ax2 = subplot(1,2,2);
box on
plot(s2_grid,x3fin,'d','LineWidth',4,'DisplayName','Numerical solution');
hold on
plot(s2_grid_exact,x3fin_true,'k-','LineWidth',2,'DisplayName','Analytic solution');
hold off
legend('FontSize',13,'Interpreter','latex');
xlabel('$s_{2}$','FontSize',15,'Interpreter','latex');
ylabel('$\mathbf{x}_{3}$','FontSize',15,'Interpreter','latex');
title(['$\mathbf{x}_{3}(t=',num2str(solution.timedep.dtime(end)),',s_{3})$'],'FontSize',14,'Interpreter','latex');
set(gca,'XLim',[min(s2_grid),max(s2_grid)]);
set(gca,'TickLabelInterpreter','latex');
ax2.Position = [0.57,0.15,0.4,0.72];
%saveas(fig1,'Ch6_ExB_1D_State','epsc');

% Plot 2D state value at final time
fig2 = figure('Position',[200 50 900 350]); 
sgtitle(['2D PDE State at Final Time, $\mathbf{x}_{4}(t=',num2str(solution.timedep.dtime(end)),',s_{1},s_{2})$'],'FontSize',15,'Interpreter','latex');
set(gcf, 'Color', 'w');
ax1 = subplot(1,2,1);
surf(s1_grid,s2_grid,x4fin','FaceAlpha',0.75,'Linestyle','--','FaceColor','interp');
h = colorbar;
%colormap jet
box on
title('Numerical solution','FontSize',14,'Interpreter','latex');
xlabel('$s_{1}$','FontSize',15,'Interpreter','latex');
ylabel('$s_{2}$','FontSize',15,'Interpreter','latex');
zlabel('$\mathbf{x}_{4}$','FontSize',15,'Interpreter','latex');
set(gca,'XLim',[min(s1_grid),max(s1_grid)]);
set(gca,'YLim',[min(s2_grid),max(s2_grid)]);
set(gca,'TickLabelInterpreter','latex');
ax1.Position = [0.08,0.11,0.3,0.72];
ax2 = subplot(1,2,2);
surf(s1_grid_exact(1:2:end),s2_grid_exact(1:2:end),x4fin_true(1:2:end,1:2:end)','FaceAlpha',0.75,'Linestyle','--','FaceColor','interp');
h = colorbar;
%colormap jet
box on
title('Analytic solution','FontSize',14,'Interpreter','latex');
xlabel('$s_{1}$','FontSize',15,'Interpreter','latex');
ylabel('$s_{2}$','FontSize',15,'Interpreter','latex');
zlabel('$\mathbf{x}_{4}$','FontSize',15,'Interpreter','latex');
set(gca,'XLim',[min(s1_grid),max(s1_grid)]);
set(gca,'YLim',[min(s2_grid),max(s2_grid)]);
set(gca,'TickLabelInterpreter','latex');
ax2.Position = [0.57,0.11,0.3,0.72];
%saveas(fig2,'Ch6_ExB_2D_State','epsc');



%% 6.4.3: DDE Simulation
clear DDE uinput; close all
% Declare the DDE to simulate
DDE.A0 = [-1 2;0 1];      DDE.Ai{1} = [.6 -.4; 0 0];
DDE.Ai{2} = [0 0; 0 -.5]; DDE.B1 = [1;1];
DDE.B2 = [0;1];           DDE.C1 = [1 0;0 1;0 0];
DDE.D12=[0;0;.1];         DDE.tau = [1,2];

% Declare initial conditions and disturbance
syms sx st sy
uinput.x = [0,0];
uinput.w(1) = -4*st-4;
uinput.u(1) = 0;
% Declare discretization and temporal integration options
opts.plot = 'yes';
opts.N = 8;
opts.tf = 1;
opts.Norder = 2;
opts.dt = 1e-3;
% Simulate and extract solution
solution = PIESIM(DDE,opts,uinput);
tval = solution.timedep.dtime;
xval = solution.timedep.ode;
zval = solution.timedep.regulated;

% Plot simulated evolution of the DDE state
fig1 = figure('Position',[200 50 900 350]); 
plot(tval,xval,'-o','LineWidth',2.5,'MarkerSize',5,'MarkerIndices',1:50:length(solution.timedep.dtime));
ax = gca;
set(ax,'XTick',[0,solution.timedep.dtime(100:100:end)]);
ax.TickLabelInterpreter = 'latex';
lgd1 = legend('$x_{1}$','$x_{2}$','Interpreter','latex','FontSize',13);
lgd1.Location = 'northeast';
title('Simulated evolution of DDE states, $x$, without state feedback control','Interpreter','latex','FontSize',14);
ylabel('$x_{1}(t)$, $x_{2}(t)$','Interpreter','latex','FontSize',14);
xlabel('$t$','FontSize',14,'Interpreter','latex');
%saveas(fig1,'Ch6_ExC_DDE_State','epsc');



%% 6.4.4 PIE Simulation
% Synthesize an optimal controller for the DDE in PIE representation
PIE = convert_PIETOOLS_DDE(DDE,'pie');
[~, K, gam] = lpiscript(PIE,'hinf-controller','light');
PIE_CL = closedLoopPIE(PIE,K);
% Declare discretization and temporal integration options
opts.plot = 'yes';
opts.N = 8;
opts.tf = 1;
opts.Norder = 2;
opts.dt = 1e-3;
% Declare initial conditions and disturbance
clear uinput;    syms st
uinput.w(1) = -4*st-4;
% Declare regularity of PIE state
ndiff = [0, 2];
% Simulate and extract solution
solution = PIESIM(PIE_CL,opts,uinput,ndiff);
tval = solution.timedep.dtime;
xval = solution.timedep.ode;
zval = solution.timedep.regulated;

% Plot simulated evolution of the DDE state
fig2 = figure('Position',[200 50 900 350]); 
plot(tval,xval,'-o','LineWidth',2.5,'MarkerSize',5,'MarkerIndices',1:50:length(solution.timedep.dtime));
ax = gca;
set(ax,'XTick',[0,solution.timedep.dtime(100:100:end)]);
ax.TickLabelInterpreter = 'latex';
lgd2 = legend('$x_{1}$','$x_{2}$','Interpreter','latex','FontSize',13);
lgd2.Location = 'northeast';
title('Simulated evolution of DDE states, $x$, with fundamental state feedback control','Interpreter','latex','FontSize',14);
ylabel('$x_{1}(t)$, $x_{2}(t)$','Interpreter','latex','FontSize',14);
xlabel('$t$','FontSize',14,'Interpreter','latex');
%saveas(fig2,'Ch6_ExD_DDE_State','epsc');
