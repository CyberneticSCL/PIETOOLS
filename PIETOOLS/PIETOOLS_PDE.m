close all; clc; clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_PDE.m     PIETOOLS 2022
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
%%
% % --- Example Library Option (See User Manual, Section 14) ---
 PDE = examples_PDE_library_PIETOOLS;
%PDE = examples_PDE_library_PIETOOLS(15,'batch');
% PDE = examples_PDE_library_PIETOOLS(15,'terms');

%% --- Manual Declaration Option --- 
% To use this example, comment lines 43
% and 45 and uncomment line 44
% pvar s t
%   A1=[0 1; 2 0]; A0=[0 0; 0 -2];
% PDE =sys();
% x1=state('pde');x2=state('pde');x=[x1;x2];
% eq_dyn=diff(x,t)==A0*x+A1*diff(x,s);
% eq_bc=[subs(x2,s,0)==0;subs(x1,s,1)==0];
% PDE=addequation(PDE,[eq_dyn;eq_bc]);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Convert to a PIE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PDE= initialize_PIETOOLS_PDE(PDE);
%PIE=convert(PDE,'pie');
PIE = convert_PIETOOLS_PDE(PDE);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Analysis or Control script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % --- Specify settings ---
%settings = settings_PIETOOLS_heavy;
settings = lpisettings('veryheavy');
settings.sos_opts.solver='sedumi';    % Solver to use
settings.eppos = 1e-4;                % Positivity of Lyapunov Function with respect to real-valued states
settings.eppos2 = 1*1e-6;             % Positivity of Lyapunov Function with respect to spatially distributed states
settings.epneg = 0;                   % Negativity of Derivative of Lyapunov Function in both ODE and PDE state -  >0 if exponential stability desired

% % --- Prompt for settings and choose executive automatically based on the example ---
%  PIETOOLS_auto_execute

% % --- Manually run desired executives ---
% [prog,P] = lpisolve(PIE,settings,'hinf-observer');
% [prog, P] = PIETOOLS_stability(PIE,settings);
% [prog, P] = PIETOOLS_stability_dual(PIE,settings);
%[prog, P, gamma] = PIETOOLS_Hinf_gain(PIE,settings);
% [prog, P, gamma] = PIETOOLS_Hinf_gain_dual(PIE,settings);
% [prog, K, gamma, P, Z] = PIETOOLS_Hinf_control(PIE,settings);
% [prog, L, gamma, P, Z] = PIETOOLS_Hinf_estimator(PIE,settings);
% [prog, Wo, gamma] = PIETOOLS_H2_norm_o(PIE,settings);
% [prog, Wc, gamma] = PIETOOLS_H2_norm_c(PIE,settings);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4: Simulation (See User Manual, Chapter 16 or xPIESIM/solver_PIESIM.m for more examples)
% Only works for PDE examples in batch input format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if PIE.dim ~=2
    opts.tf=10;
    opts.plot='yes';
    uinput.ic.ODE=1;

    [solution, grid] = PIESIM(PDE,opts,uinput);

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
        X = repmat(grid.phys(:,1),1,length(t));
        n1=0;
        if PIE.dim ~=2
        pde_sol = solution.timedep.pde;
        else
            if ~isempty(solution.timedep.pde{1})
        n=size(solution.timedep.pde{1});
        n1=n(2);
        pde_sol = solution.timedep.pde{1};
            end
            if (length(solution.timedep.pde)==2)
            if ~isempty(solution.timedep.pde{2})
        k=floor(length(grid.phys)/2)+1;
        n=size(solution.timedep.pde{2});
        pde_sol(1:n(1),n1+1:n1+n(3),1:n(4)) = solution.timedep.pde{2}(1:n(1),k,1:n(3),1:n(4));
            end
            end
        end
        t = repmat(t,size(pde_sol,1),1);
        figure;
        for i=1:size(pde_sol,2)
            Z = squeeze(pde_sol(:,i,:));
            ax(i)=subplot(size(pde_sol,2),1,i);
            surf(t,grid.phys(:,1),Z,'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
            h=colorbar ;
            colormap jet
            box on
            ylabel(h,'$|\mathbf{x}(t,s)|$','interpreter', 'latex','FontSize',15)
            set(gcf, 'Color', 'w');
            xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$s$','FontSize',15,'Interpreter','latex');
            zlabel('$\mathbf{x}(t,s)$','FontSize',15,'Interpreter','latex');
        end
        subplot(ax(1));
        title('Time evolution of open-loop PDE states, x, plotted against space, s');
    end
