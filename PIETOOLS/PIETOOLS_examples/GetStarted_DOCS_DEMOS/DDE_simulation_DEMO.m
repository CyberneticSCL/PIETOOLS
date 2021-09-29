%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DDE_simulation_DEMO.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This is a demo file for simulation of DDEs using PIESIM code.
% It perfoms the following functions:
% 1) Sets up the DDE problem (either by calling an example from the library
% or setting it up manually)
% 2) Finds a PIE (constructs closed-loop observer or controller if
% specified)
% 2) Call the executive solver for PIESIM (exectuive_PIESIM.m)

% The executive rotuine performs the following functions:
% 1) Checks if all inputs are defined and (if not) setting up default
% options
% 1) Converts the DDE to PIE problem
% 2) Discretizes the PIE operators with Chebyshev polynomials
% 3) Temporally integrates the corresponding spatially-discretized
% equations (3 temporal schemes are possible: BDF, Gauss integration, and
% analytical integraiton - see options in solver_PIESIM.m)
% 4) Outputs and plots the PDE and ODE solutions at a final time and

clear;
clc;
close all;
format long;
pvar s theta;
Hinf_control=0;

%--------------------------------------------------------------
% Simulation Examples:Uncomment one of the following 3 examples to run the simulations
%--------------------------------------------------------------
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Example B.2 from [4], adapted from [7]
% % % % \dot x=.1x(t)-x(t-tau)+w(t)+u(t)
% % % % % z(t)=x(t)+w(t)
% Hinf_control=1
% DDE.A0=[2 1;0 -1];%
% DDE.Ai{1}=[-1 0; -1 1];%
% DDE.B1=[-.5;1]; DDE.B2=[3;1];
% DDE.C1=[1 -.5;0 0];
% DDE.D12=[0;1];
% DDE.tau(1) = .3;
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Example B.3 from [4], adapted from [6]
Hinf_control=1
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
% Define simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the Chebyshev polynomial discretization order 
opts.N=8;

% final simulation time
opts.tf=30;

% Integration scheme
% opts.intScheme = 1 - Backward Difference Formula (BDF)
opts.intScheme=1;
opts.Norder=2; % order of truncation term in BDF scheme
opts.dt=0.01; % gap between time steps
opts.Nsteps=floor(opts.tf/opts.dt); % number of time steps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert DDE to PIE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DDE=initialize_PIETOOLS_DDE(DDE); 
DDF=minimize_PIETOOLS_DDE2DDF(DDE);
PIE=convert_PIETOOLS_DDF2PIE(DDF);

% extract relevant PIE dimensions
psize.nu=PIE.Tu.dim(1,2);
psize.nw=PIE.Tw.dim(1,2);
psize.nx=PIE.T.dim(1,1);
psize.nf=0;
ns=PIE.T.dim(2,1);
psize.N=opts.N;
psize.n=ns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define inputs and IC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms sx st;
uinput.ifexact=false;
uinput.a=-1;
uinput.b=0;
if (psize.nw>0)
    if ~isfield(uinput,'w')
        uinput.w(1:psize.nw)=sin(st);
        uinput.w(1:psize.nw)=1+0*st;
        uinput.w(1:psize.nw)=sin(st);
        uinput.w(1:psize.nw)=sinc(st);
    end
end
n_diff = PIE.T.dim(2,1);
uinput.ic.ODE(1:psize.nx)=1;
uinput.ic.PDE(1:ns)=1;

plot_no_control = 1; % this setting first runs simulation without controller for same IC

if plot_no_control % run simulation without controller
    solutionA=executive_PIESIM(DDE,opts,uinput);
    close all;
    
end

if Hinf_control
    % find the controller and find closed Loop PIE
    settings_PIETOOLS_custom
    settings.sos_opts.solver = 'sedumi';
    settings.eppos = 1e-4;      % Positivity of Lyapunov Function with respect to real-valued states
    settings.eppos2 = 0*1e-4;   % Positivity of Lyapunov Function with respect to spatially distributed states
    settings.epneg = 0*1e-5;    % Negativity of Derivative of Lyapunov Function
    [prog, K, gamma, P, Z] = executive_PIETOOLS_Hinf_control(PIE,settings);
    PIE = closedLoopPIE(PIE,K);
    
    solutionB=executive_PIESIM(PIE,opts,uinput,n_diff);
    close all;
end

for i=1:psize.nx
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


