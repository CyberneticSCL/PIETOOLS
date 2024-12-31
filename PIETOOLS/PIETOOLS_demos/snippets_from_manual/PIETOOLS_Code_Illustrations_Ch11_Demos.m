% This document provides the codes used to plot the simulation results from
% the various demos, as they appear in Chapter 11 the manual. We refer to 
% the demo codes and Chapter 11 of the manual for more information.

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
% DJ, 12/30/2024: Initial coding;
clear

%% 11.1 DEMO 1: Stability and Simulation
% Run the Demo
DEMO1_stability_and_simulation

% % Plot the evolution of the two PDE state variables
plot_indcs = floor(linspace(1,size(phi1,2),100)); 
tplot = tval(plot_indcs);           % Only plot at select times
colors = {'b','g','m','r','k','c','r','y','o'};     % colors for the plot
fig1 = figure('Position',[200 150 1000 400]); 
set(gcf, 'Color', 'w');
subplot(1,2,1);
box on
surf(tplot,grids.phys,phi1(:,plot_indcs),'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
h = colorbar;
colormap jet    
subplot(1,2,2);
box on
surf(tplot,grids.phys,phi2(:,plot_indcs),'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
h = colorbar;
colormap jet
% Clean up the figure
sgtitle('PDE state evolution','Interpreter','latex','FontSize',16)
ax1 = subplot(1,2,1);       ax1.TickLabelInterpreter = 'latex';
ax1.Position = [0.06,0.1,0.33,0.75];
box on
title('$\phi_{1}(t,s)=\partial_{s} \mathbf{x}(t,s)$','Interpreter','latex','FontSize',14);
xlabel('$t$','FontSize',15,'Interpreter','latex');  
ylabel('$s$','FontSize',15,'Interpreter','latex');
zlabel('$\phi$','FontSize',15,'Interpreter','latex');
ax2 = subplot(1,2,2);       ax2.TickLabelInterpreter = 'latex';
ax2.Position = [0.57,0.1,0.33,0.75];
box on
title('$\phi_{2}(t,s)=\partial_{t} \mathbf{x}(t,s)$','Interpreter','latex','FontSize',14);
xlabel('$t$','FontSize',15,'Interpreter','latex');  
ylabel('$s$','FontSize',15,'Interpreter','latex');
zlabel('$\phi$','FontSize',15,'Interpreter','latex');
%saveas(fig1,'DEMO1_StateEvolution_Plot','epsc');

% % Plot the ODE state and regulated output versus time
fig2 = figure('Position',[200 150 1000 350]);
set(gcf, 'Color', 'w');
subplot(1,2,1)
box on
plot([0,tplot],[0.5,xval(1,plot_indcs)],[colors{2},'-'],'LineWidth',1.5,'DisplayName','$x(t)$');
subplot(1,2,2)
box on
hold on
plot([0,tplot],[1,zval(1,plot_indcs)],[colors{1},'-'],'LineWidth',1.5,'DisplayName','$r(t)$');
plot([0,tplot],[0,wval(plot_indcs)],[colors{4},'-.'],'LineWidth',1.0,'DisplayName','$w(t)$');
hold off
% Clean up the figure
sgtitle('ODE state and regulated output evolution','Interpreter','latex','FontSize',15)
ax1 = subplot(1,2,1);   ax1.XLim = [0,opts.tf];
ax1.TickLabelInterpreter = 'latex';
title('ODE state $x(t)$','Interpreter','latex','FontSize',13)
xlabel('$t$','FontSize',15,'Interpreter','latex');
ylabel('$x$','FontSize',15,'Interpreter','latex');
ax2 = subplot(1,2,2);   ax2.XLim = [0,opts.tf];
ax2.TickLabelInterpreter = 'latex';
title('Regulated output $r(t)$ and disturbance $w(t)$','Interpreter','latex','FontSize',13)
legend('Interpreter','latex','Location','northeast','FontSize',10.5);
xlabel('$t$','FontSize',15,'Interpreter','latex');    
ylabel('$r$','FontSize',15,'Interpreter','latex');
%saveas(fig2,'DEMO1_OutputEvolution_Plot','epsc');






%% 11.5 DEMO 5: Hinf Optimal Estimator
% Run the Demo
DEMO5_Hinf_optimal_estimator

% % Plot the simulated true and estimated state, as well as the error in
% % the state estimate at several grid points

% Set options for the plot
grid_idcs = [1,5,7];        % Only plot at first, fifth and seventh grid points
plot_indcs = floor(logspace(0,log(0.5*opts.tf/opts.dt)/log(10),40)); 
%plot_indcs = floor(linspace(1,opts.tf/opts.dt,66)); 
tplot = tval(plot_indcs);           % Only plot at select times
colors = {'b','g','m','r','k','c','r','y','o'};     % Colors for the plot

% Plot evolution of actual and estimated state at grid points
fig1 = figure();
set(gcf, 'Color', 'w');
box on
for j=grid_idcs
    s_pos = num2str(grid.phys(j));  % Position associated to grid index.
    subplot(1,2,1)
    box on
    hold on
    plot(tplot,x_act(j,plot_indcs),[colors{j},'-'],'LineWidth',2,'DisplayName',['$\mathbf{x}(s=',s_pos,')$']);
    plot(tplot,x_est(j,plot_indcs),[colors{j},'o--'],'LineWidth',1.5,'DisplayName',['$\mathbf{\hat{x}}(s=',s_pos,')$']);
    hold off
    
    subplot(1,2,2)
    box on
    hold on
    loglog(tplot,(x_act(j,plot_indcs)-x_est(j,plot_indcs)),[colors{j},'o-'],'LineWidth',1.5,'DisplayName',['$\mathbf{e}(s=',s_pos,')$']);
    hold off
end

% Clean up the figure
ax1 = subplot(1,2,1);   ax1.XScale = 'log';     
ax1.XTickLabels = {'0.001';'0.01';'0.1';'1'};
ax1.TickLabelInterpreter = 'latex';
lgd1 = legend('Interpreter','latex');                lgd1.FontSize = 10.5;
lgd1.Location = 'northwest';
xlabel('$t$','FontSize',15,'Interpreter','latex');  ylabel('$\mathbf{x}$','FontSize',15,'Interpreter','latex');
title('PDE state value $\mathbf{x}(t)$ and estimate $\mathbf{\hat{x}}(t)$','Interpreter','latex','FontSize',15);
ax1 = subplot(1,2,2);   ax1.XScale = 'log';     
ax1.XTickLabels = {'0.001';'0.01';'0.1';'1'};
ax1.TickLabelInterpreter = 'latex';
lgd2 = legend('Interpreter','latex');                lgd2.FontSize = 10.5;
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$\mathbf{e}$','FontSize',15,'Interpreter','latex');
title('Error $\mathbf{e}=\mathbf{\hat{x}}(t)-\mathbf{x}(t)$','Interpreter','latex','FontSize',15);
fig1.Position = [200 150 1000 450];

%saveas(fig1,'DEMO5_EstimateError_Plot','epsc');



%% 11.6 DEMO 6: Hinf Optimal Control
% Run the Demo
DEMO6_Hinf_optimal_control

% % Compute the initial value of the state at each grid point.
% (check that d^2/ds^2 x0_sym = sin(pi/2 *s);
x0_sym = -(4/sym(pi)^2)*sin((sym(pi)/2)*sx);

% % Plot the open- and closed-loop state at different times
plot_indcs = [1,5,10,50,100];
colors = {'b','g','m','r','k','c','r','y','o'};     % colors for the plot
figure('Position',[200 150 1000 450]); 
set(gcf, 'Color', 'w');
box on
% Plot initial state value
subplot(1,2,1);
plot(grid.phys,subs(x0_sym,sx,grid.phys),[colors{1},'-'],'LineWidth',2,'DisplayName',['$\mathbf{x}_{ol}(t=0)$']);
subplot(1,2,2);
plot(grid.phys,subs(x0_sym,sx,grid.phys),[colors{1},'-'],'LineWidth',2,'DisplayName',['$\mathbf{x}_{cl}(t=0)$']);
% Plot state at later times
for j=1:length(plot_indcs)
    plot_indx = plot_indcs(j);
    t_pos = num2str(tval(plot_indx));
    
    subplot(1,2,1);
    hold on
    plot(grid.phys,x_OL(:,plot_indx),[colors{j+1},'-'],'LineWidth',2,'DisplayName',['$\mathbf{x}_{ol}(t=',t_pos,')$']);
    hold off
    
    subplot(1,2,2);
    hold on
    plot(grid.phys,x_CL(:,plot_indx),[colors{j+1},'-'],'LineWidth',2,'DisplayName',['$\mathbf{x}_{cl}(t=',t_pos,')$']);
    hold off
end
% Clean up the figure
sgtitle('PDE state $\mathbf{x}(t,s)$ at various times','Interpreter','latex','FontSize',15)
ax1 = subplot(1,2,1);       ax1.TickLabelInterpreter = 'latex';
box on
title('Open-loop','Interpreter','latex','FontSize',13);
legend('Interpreter','latex','Location','northwest','FontSize',10.5,'Direction','reverse');
xlabel('$s$','FontSize',15,'Interpreter','latex');  ylabel('$\mathbf{x}_{ol}$','FontSize',15,'Interpreter','latex');
ax2 = subplot(1,2,2);       ax2.TickLabelInterpreter = 'latex';
box on
title('Closed-loop','Interpreter','latex','FontSize',13);
legend('Interpreter','latex','Location','northwest','FontSize',10.5,'Direction','reverse');
xlabel('$s$','FontSize',15,'Interpreter','latex');    ylabel('$\mathbf{x}_{cl}$','FontSize',15,'Interpreter','latex');

% % Plot the regulated output versus time
plot_indcs = floor(logspace(0,log(opts.tf/opts.dt)/log(10),50));
tplot = tval(plot_indcs);
figure('Position',[200 150 1000 450]);
set(gcf, 'Color', 'w');
box on
subplot(1,2,1)
box on
hold on
plot(tplot,z_OL(1,plot_indcs),'-','LineWidth',1.5,'DisplayName','$z_{ol}(t)$');
plot(tplot,w(plot_indcs),'-.','LineWidth',1.0,'DisplayName','$w(t)$');
hold off
subplot(1,2,2)
box on
hold on
plot(tplot,z_CL(1,plot_indcs),'-','LineWidth',1.5,'DisplayName','$z_{cl}(t)$');
plot(tplot,w(plot_indcs),'-.','LineWidth',1.0,'DisplayName','$w(t)$');
hold off
% Clean up the figure
sgtitle('Regulated output $z(t)=\int_{0}^{1}\mathbf{x}(t,s)ds+w(t)$ and disturbance $w(t)$','Interpreter','latex','FontSize',15)
ax1 = subplot(1,2,1);   ax1.XScale = 'log';  ax1.XLim = [opts.dt,opts.tf];
ax1.TickLabelInterpreter = 'latex';
title('Open-loop','Interpreter','latex','FontSize',13)
legend('Interpreter','latex','Location','northeast','FontSize',10.5);
xlabel('$t$','FontSize',15,'Interpreter','latex');
ylabel('$z_{ol}$','FontSize',15,'Interpreter','latex');
ax2 = subplot(1,2,2);   ax2.XScale = 'log';  ax2.XLim = [opts.dt,opts.tf];
ax2.TickLabelInterpreter = 'latex';
title('Closed-loop','Interpreter','latex','FontSize',13)
legend('Interpreter','latex','Location','northeast','FontSize',10.5);
xlabel('$t$','FontSize',15,'Interpreter','latex');    
ylabel('$z_{cl}$','FontSize',15,'Interpreter','latex');

% % Plot the control effort
plot_indcs = floor(logspace(0,log(opts.tf/opts.dt)/log(10),50));
tplot = tval(plot_indcs);
figure('Position',[200 150 1000 450]);
set(gcf, 'Color', 'w');
box on
plot(tplot,u_CL(1,plot_indcs),[colors{7},'-'],'LineWidth',1.5,'DisplayName','$u_{cl}(t)$');
title('Closed-loop control effort $u_{cl}(t)=\mathcal{K}\mathbf{x}_{f,cl}(t)$','Interpreter','latex','FontSize',15)
ax3 = gca;   ax3.XScale = 'log';  ax3.XLim = [opts.dt,opts.tf];
ax3.TickLabelInterpreter = 'latex';
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$u_{cl}$','FontSize',15,'Interpreter','latex');

%saveas(figure(1),'DEMO6_State_Plot','epsc');
%saveas(figure(2),'DEMO6_Output_Plot','epsc');
%saveas(figure(3),'DEMO6_Control_Plot','epsc');



%% 11.7 DEMO 7: Observer Based Control
DEMO7_observer_based_control

% % Plot the open- and closed-loop state evolution
plot_indcs = floor(logspace(0,log(0.5*opts.tf/opts.dt)/log(10),40)); 
tplot = tval(plot_indcs);           % Only plot at select times
colors = {'b','g','m','r','k','c','r','y','o'};     % colors for the plot
figure('Position',[200 150 1000 400]); 
set(gcf, 'Color', 'w');
subplot(1,2,1);
box on
surf(tplot,grid.phys,x_OL(:,plot_indcs),'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
h = colorbar;
colormap jet    
subplot(1,2,2);
box on
surf(tplot,grid.phys,x_CL(:,plot_indcs),'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
h = colorbar;
colormap jet
% Clean up the figure
sgtitle('PDE state evolution $\mathbf{x}(t,s)$','Interpreter','latex','FontSize',16)
ax1 = subplot(1,2,1);       ax1.TickLabelInterpreter = 'latex';
ax1.XScale = 'log';
ax1.Position = [0.06,0.1,0.33,0.75];
box on
title('Open-loop','Interpreter','latex','FontSize',14);
xlabel('$t$','FontSize',15,'Interpreter','latex');  
ylabel('$s$','FontSize',15,'Interpreter','latex');
zlabel('$\mathbf{x}_{ol}$','FontSize',15,'Interpreter','latex');
ax2 = subplot(1,2,2);       ax2.TickLabelInterpreter = 'latex';
ax2.XScale = 'log';
ax2.Position = [0.57,0.1,0.33,0.75];
box on
title('Closed-loop','Interpreter','latex','FontSize',14);
xlabel('$t$','FontSize',15,'Interpreter','latex');  
ylabel('$s$','FontSize',15,'Interpreter','latex');
zlabel('$\mathbf{x}_{cl}$','FontSize',15,'Interpreter','latex');

% % Plot the regulated output versus time
plot_indcs = floor(logspace(0,log(opts.tf/opts.dt)/log(10),50));
tplot = tval(plot_indcs);
figure('Position',[200 150 1000 400]);
set(gcf, 'Color', 'w');
subplot(1,2,1)
box on
hold on
plot(tplot,z_OL(1,plot_indcs),'-','LineWidth',1.5,'DisplayName','$z_{ol}(t)$');
plot(tplot,w(plot_indcs),'-.','LineWidth',1.0,'DisplayName','$w(t)$');
hold off
subplot(1,2,2)
box on
hold on
plot(tplot,z_CL(1,plot_indcs),'-','LineWidth',1.5,'DisplayName','$z_{cl}(t)$');
plot(tplot,w(plot_indcs),'-.','LineWidth',1.0,'DisplayName','$w(t)$');
hold off
% Clean up the figure
sgtitle('Regulated output $z(t)=\int_{0}^{1}\mathbf{x}(t,s)ds$ and disturbance $w(t)$','Interpreter','latex','FontSize',15)
ax1 = subplot(1,2,1);   ax1.XScale = 'log';  ax1.XLim = [opts.dt,opts.tf];
ax1.TickLabelInterpreter = 'latex';
title('Open-loop','Interpreter','latex','FontSize',13)
legend('Interpreter','latex','Location','northwest','FontSize',10.5);
xlabel('$t$','FontSize',15,'Interpreter','latex');
ylabel('$z_{ol}$','FontSize',15,'Interpreter','latex');
ax2 = subplot(1,2,2);   ax2.XScale = 'log';  ax2.XLim = [opts.dt,opts.tf];
ax2.TickLabelInterpreter = 'latex';
title('Closed-loop','Interpreter','latex','FontSize',13)
legend('Interpreter','latex','Location','northeast','FontSize',10.5);
xlabel('$t$','FontSize',15,'Interpreter','latex');    
ylabel('$z_{cl}$','FontSize',15,'Interpreter','latex');

% % Plot the control effort
plot_indcs = floor(logspace(0,log(opts.tf/opts.dt)/log(10),50));
tplot = tval(plot_indcs);
figure('Position',[200 150 1000 350]);
set(gcf, 'Color', 'w');
box on
plot(tplot,u_CL(1,plot_indcs),[colors{7},'-'],'LineWidth',1.5,'DisplayName','$u(t)$');
title('Closed-loop control effort $u(t)=\mathcal{K}\hat{\mathbf{x}}_{f}(t)$','Interpreter','latex','FontSize',15)
ax3 = gca;   ax3.XScale = 'log';  ax3.XLim = [opts.dt,opts.tf];
ax3.TickLabelInterpreter = 'latex';
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$u$','FontSize',15,'Interpreter','latex');

%saveas(figure(1),'DEMO7_State_Plot','epsc');
%saveas(figure(2),'DEMO7_Output_Plot','epsc');
%saveas(figure(3),'DEMO7_Control_Plot','epsc');



%% 11.9 DEMO 9: L2-Gain Analysis and Simulation of 2D PDEs
% % Plot the regulated output versus time
%plot_indcs = floor(logspace(0,log(opts.tf/opts.dt)/log(10),50));
plot_indcs = round(linspace(1,size(z,2),100));
tplot = tval(plot_indcs);
figure('Position',[200 150 1000 350]);
set(gcf, 'Color', 'w');
subplot(1,2,1)
box on
plot(tplot,w(plot_indcs),'r-.','LineWidth',1.5,'DisplayName','$w(t)$');
subplot(1,2,2)
box on
plot(tplot,z(1,plot_indcs),'b-','LineWidth',1.5,'DisplayName','$z(t)$');
% Clean up the figure
sgtitle('Regulated output response to bounded disturbance','Interpreter','latex','FontSize',15)
ax1 = subplot(1,2,1);   ax1.XScale = 'linear';  ax1.XLim = [opts.dt,opts.tf];
ax1.TickLabelInterpreter = 'latex';
title('Disturbance $w(t)=20\sin(\pi t)e^{-t/2}$','Interpreter','latex','FontSize',13)
xlabel('$t$','FontSize',15,'Interpreter','latex');
ylabel('$w$','FontSize',15,'Interpreter','latex');
ax2 = subplot(1,2,2);
ax2.TickLabelInterpreter = 'latex';
title('Regulated output $z(t)=\int_{0}^{1}\int_{0}^{1}\mathbf{x}(t,s)ds$','Interpreter','latex','FontSize',13)
xlabel('$t$','FontSize',15,'Interpreter','latex');
ylabel('$z$','FontSize',15,'Interpreter','latex');

%saveas(figure(3),'DEMO9_OutputResponse_Plot','epsc');