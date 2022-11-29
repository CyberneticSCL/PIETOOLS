%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIE_simulation_SNIPPET.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code Snippet illustrating usage of PIESIM for simulation of PIEs.
% See Section 6.3.3 in Manual for a description.
% See also "DDE_simulation_SNIPPET.m".

%% Defining DDE and finding an optimal controller
DDE.A0=[-1 2;0 1]; DDE.Ai{1}=[.6 -.4; 0 0];
DDE.Ai{2}=[0 0; 0 -.5]; DDE.B1=[1;1];
DDE.B2=[0;1]; DDE.C1=[1 0;0 1;0 0];
DDE.D12=[0;0;.1]; DDE.tau=[1,2];
DDE=initialize_PIETOOLS_DDE(DDE);
DDF=minimize_PIETOOLS_DDE2DDF(DDE);
PIE=convert_PIETOOLS_DDF(DDF);
sett = lpisettings('heavy');
[prog, K, gamma, P, Z] = PIETOOLS_Hinf_control(PIE,sett);
%% constructing closed loop system
PIE = closedLoopPIE(PIE,K);
ndiff = [0, PIE.T.dim(2,1)];

%% Setting PIESIM simulation parameters
syms st;
uinput.w = -4*st-4;
uinput.u = 0;
uinput.ic.ODE = [1;0];
opts.plot='no';
opts.N=8;
opts.tf=1;
opts.intScheme=1;
opts.Norder = 2;
opts.dt=1e-3;

%% Simulating and plotting open loop system
solution_ol = PIESIM(DDE,opts,uinput);
%%
plot(solution_ol.timedep.dtime,solution_ol.timedep.ode,'--o','MarkerIndices',1:50:length(solution_ol.timedep.dtime));
ax = gca;
set(ax,'XTick',solution_ol.timedep.dtime(1:150:end));
lgd1 = legend('$x_1$','$x_2$','Interpreter','latex'); lgd1.FontSize = 10.5; 
lgd1.Location = 'northeast';
title('Time evolution of the Delay system states, x, without state feedback control');
ylabel('$x_1(t), ~~~x_2(t)$','Interpreter','latex','FontSize',15);
xlabel('t','FontSize',15,'Interpreter','latex');
%% Simulating and plotting closed loop system
opts.tf=10;
uinput.ic.PDE = [0;0]; uinput.ic.ODE = [0,0];
solution_cl=PIESIM(PIE,opts,uinput,ndiff);

plot(solution_cl.timedep.dtime,solution_cl.timedep.ode,'--o','MarkerIndices',1:50:length(solution_cl.timedep.dtime));
ax = gca;
set(ax,'XTick',solution_cl.timedep.dtime(1:150:end));
lgd1 = legend('$x_1$','$x_2$','Interpreter','latex'); lgd1.FontSize = 10.5; 
lgd1.Location = 'northeast';
title('Time evolution of the Delay system states, x, with state feedback control');
ylabel('$x_1(t), ~~~x_2(t)$','Interpreter','latex','FontSize',15);
xlabel('t','FontSize',15,'Interpreter','latex');