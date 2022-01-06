close all; clc; clear;
pvar s theta; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_PDE.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting Started: 
%
% Step 1: Declare a PDE model (Batch Format, Example file or GUI)
% Step 2: Convert to a PIE
% Step 3: Run an analysis or control script
% Step 4: Run PIESIM

% This program determines stability, hinf-norm and designs hinf-optimal
% observer of a linear coupled ODE-PDE which is defined in the format given below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Declare a PDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % --- GUI option (See User Manual, Section 12) ---
% PIETOOLS_PDE_GUI

% % --- Example Library Option (See User Manual, Section 14) ---
 PDE = examples_PDE_library_PIETOOLS;
% PDE = examples_PDE_library_PIETOOLS(5,'batch');
% PDE = examples_PDE_library_PIETOOLS(5,'terms');

% % --- Batch or Terms Declaration Option (See User Manual, Section 13) ---
% PDE_b.n0 = 0;   PDE_b.n1 = 1;   PDE_b.n2 = 0;   % state dimensions
% PDE_b.dom = [0,1];                              % spatial domain
% PDE_b.A1= 1;                                    % PDE dx/dt = A*dx/ds
% PDE_b.B = [eye(PDE_b.n1) zeros(PDE_b.n1)];      % BC x(0) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Convert to a PIE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 PIE = convert_PIETOOLS_PDE(PDE);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Analysis or Control script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % --- Specify settings ---
% settings_PIETOOLS_heavy;
 settings_PIETOOLS_light;
% settings_PIETOOLS_stripped;
% settings.sos_opts.solver='sedumi';    % Solver to use
% settings.eppos = 1e-4;                % Positivity of Lyapunov Function with respect to real-valued states
% settings.eppos2 = 1*1e-6;             % Positivity of Lyapunov Function with respect to spatially distributed states
% settings.epneg = 0;                   % Negativity of Derivative of Lyapunov Function in both ODE and PDE state -  >0 if exponential stability desired

% % --- Prompt for settings and choose executive automatically based on the example ---
 PIETOOLS_auto_execute

% % --- Manually run desired executives ---
% [prog, P] = PIETOOLS_stability(PIE,settings);
% [prog, P] = PIETOOLS_stability_dual(PIE,settings);
% [prog, P, gamma] = PIETOOLS_Hinf_gain(PIE,settings);
% [prog, P, gamma] = PIETOOLS_Hinf_gain_dual(PIE,settings);
% [prog, K, gamma, P, Z] = PIETOOLS_Hinf_control(PIE,settings);
% [prog, L, gamma, P, Z] = PIETOOLS_Hinf_estimator(PIE,settings);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Simulation (See User Manual, Chapter 16 or xPIESIM/solver_PIESIM.m for more examples)
% Only works for PDE examples in batch input format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[solution, grid] = PIESIM(PDE);

% Note: you can also specify time stepping options and user inputs, such as initial
% conditions and non-zero boundary inputs, via
% solution = executive_PIESIM(PDE, opts), 
% solution = executive_PIESIM(PDE, [], uinput), or 
% solution = executive_PIESIM(PDE, opts, uinput), please consult the user's manual

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4b: Plotting solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(ismember(fieldnames(solution),'timedep'))
    t = solution.timedep.dtime;
    X = repmat(grid.phys,1,length(t));
    pde_sol = solution.timedep.pde;
    t = repmat(t,size(pde_sol,1),1);
    figure;
    for i=1:size(pde_sol,2)
        Z = squeeze(pde_sol(:,i,:));
        ax(i)=subplot(size(pde_sol,2),1,i);
        surf(t, X, Z);
        shading interp;
        xlabel('t');
        ylabel('s');
        zlabel(['x_',num2str(i),'(t,s)'],'Interpreter','tex');
    end
    subplot(ax(1));
    title('Time evolution of PDE states, x, plotted against space, s');
end