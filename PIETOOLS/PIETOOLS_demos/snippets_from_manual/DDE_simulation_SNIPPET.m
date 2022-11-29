%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DDE_simulation_SNIPPET.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code Snippet illustrating usage of PIESIM for simulation of DDEs.
% See Section 6.3.2 in Manual for a description.

% This code perfoms the following functions:
% 1) Sets up the DDE problem (either by calling an example from the library
% or setting it up manually)
% 2) Finds a PIE (constructs closed-loop observer or controller if
% specified)
% 2) Call the executive solver for PIESIM (exectuive_PIESIM.m)

% The executive routine performs the following functions:
% 1) Checks if all inputs are defined and (if not) setting up default
% options
% 1) Converts the DDE to PIE problem
% 2) Discretizes the PIE operators with Chebyshev polynomials
% 3) Temporally integrates the corresponding spatially-discretized
% equations (3 temporal schemes are possible: BDF, Gauss integration, and
% analytical integraiton - see options in solver_PIESIM.m)
% 4) Outputs and plots the PDE and ODE solutions at a final time and

clear; clc;
close all;
format long;
pvar s theta;
plot_no_control = 1; % this setting is 1, first run simulation without controller for same IC
%--------------------------------------------------------------
% Simulation Examples:Uncomment one of the following 3 examples to run the simulations
%--------------------------------------------------------------
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Example B.2 from [4], adapted from [7]
% % % % \dot x=.1x(t)-x(t-tau)+w(t)+u(t)
% % % % % z(t)=x(t)+w(t)
% DDE.A0=[2 1;0 -1];%
% DDE.Ai{1}=[-1 0; -1 1];%
% DDE.B1=[-.5;1]; DDE.B2=[3;1];
% DDE.C1=[1 -.5;0 0];
% DDE.D12=[0;1];
% DDE.tau(1) = .3;
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Example B.3 from [4], adapted from [6]
   DDE.A0=[-1 2;0 1];%
   DDE.Ai{1}=[.6 -.4; 0 0];%
   DDE.Ai{2}=[0 0; 0 -.5];%
   DDE.B1=[1;1]; DDE.B2=[0;1];
   DDE.C1=[1 0;0 1;0 0];
   DDE.D12=[0;0;.1];
   DDE.tau(1) = 1;DDE.tau(2)=2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % Wind Tunnel Problem from [4], adapted from [5]
% Hinf_control=1
% aM=1/1.964; kM=-.0117;zM=.8;wM=6;
% DDE.A0=[-aM 0 0; 0 0 1; 0 -wM^2 -2*zM*wM];
%   Ad{1}=[0 kM*aM 0; 0 0 0; 0 0 0];
%  DDE.C0=[1 0 0; 0 1 0;0 0 0]; DDE.C{1}=[0 0 0;0 0 0;0 0 0];
%  DDE.B2=[0;0;wM^2]; %B in paper
%  DDE.B1=[1 0;0 0;0 10]; % D in paper
%  DDE.tau(1)=.33;
%  DDE.D1=[0 0;0 0;0 0];
%  DDE.D2=[0;0;.1];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert DDE to PIE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DDE = initialize_PIETOOLS_DDE(DDE); 
DDF = convert_PIETOOLS_DDE(DDE,'ddf');
PIE = convert_PIETOOLS_DDF(DDF,'pie');

if plot_no_control % run simulation without controller
    solutionA = PIESIM(DDE); 
end

if Hinf_control
    % use default settings with low order polynomials for controller
    settings = lpisettings('light');
    settings.options1.sep = 1; % this is necessary for invertible P operator
    
    % define LMI solver parameters
    settings.sos_opts.solver = 'sedumi';
    settings.eppos = 1e-4;      % Positivity of Lyapunov Function with respect to real-valued states
    settings.eppos2 = 0;   % Positivity of Lyapunov Function with respect to spatially distributed states
    settings.epneg = 0;    % Negativity of Derivative of Lyapunov Function
    
    [prog, K, gamma, P, Z] = PIETOOLS_Hinf_control(PIE,settings);
    PIE = closedLoopPIE(PIE,K); % find closed loop system with controller
    ndiff = [0, PIE.T.dim(2,1)]; % this means all distributed states are differentiable 
    solutionB = PIESIM(PIE,ndiff);
end
close all;


for i=1:PIE.T.dim(1,1)
    labels{i} = ['x_',num2str(i)];
end


figure;
plot(solutionA.timedep.dtime, solutionA.timedep.ode,'-x'); hold on;
title('Time evolution of DDE states - without controller feedback');
xlabel('t');
ylabel('x_i(t)');
legend(labels);

figure;
plot(solutionB.timedep.dtime, solutionB.timedep.ode,'-x'); hold on;
title('Time evolution of DDE states - with controller feedback');
xlabel('t');
ylabel('x_i(t)');
legend(labels);


