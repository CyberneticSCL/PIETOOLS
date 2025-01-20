close all; clc; clear; clear stateNameGenerator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_PDE.m     PIETOOLS 2024
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
% The example library contains 52 examples of 1D and 2D ODE-PDE systems,
% each with an associated default executive to run. Call the library with
% an index 1 through 52 as argument to extract that example, and run the
% associated executive, if desired.
PDE = examples_PDE_library_PIETOOLS(1);

%% --- Manual Declaration Option --- 
% % To use this example, comment line 28, and
% % uncomment line 63
% pvar s t
% A1 = [0 1; 2 0];        A0 = [0 0; 0 -2];
% x1 = pde_var(s,[0,1]);  x2 = pde_var(s,[0,1]);    x = [x1;x2];
% eq_dyn = diff(x,t)==A0*x+A1*diff(x,s);
% eq_bc = [subs(x2,s,0)==0; subs(x1,s,1)==0];
% PDE = initialize([eq_dyn;eq_bc]);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Convert to a PIE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PIE=convert(PDE,'pie');
PIE = convert_PIETOOLS_PDE(PDE);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Analysis or Control script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % --- Specify settings ---
%settings = settings_PIETOOLS_heavy;
settings = lpisettings('heavy');
settings.sos_opts.solver='sedumi';    % Solver to use
settings.eppos = 1e-4;                % Positivity of Lyapunov Function with respect to real-valued states
settings.eppos2 = 1*1e-6;             % Positivity of Lyapunov Function with respect to spatially distributed states
settings.epneg = 0;                   % Negativity of Derivative of Lyapunov Function in both ODE and PDE state -  >0 if exponential stability desired

% % --- Prompt for settings and choose executive automatically based on the example ---
% PIETOOLS_auto_execute

% % --- Manually run desired executives ---
% [prog,P] = lpisolve(PIE,settings,'hinf-observer');
% [prog, P] = PIETOOLS_stability(PIE,settings);
% [prog, P] = PIETOOLS_stability_dual(PIE,settings);
%[prog, P, gamma] = PIETOOLS_Hinf_gain(PIE,settings);
% [prog, P, gamma] = PIETOOLS_Hinf_gain_dual(PIE,settings);
% [prog, K_control, gamma, P, Z] = PIETOOLS_Hinf_control(PIE,settings);
% [prog, L_estimator, gamma, P, Z] = PIETOOLS_Hinf_estimator(PIE,settings);
% [prog, Wo, gamma] = PIETOOLS_H2_norm_o(PIE,settings);
% [prog, Wc, gamma] = PIETOOLS_H2_norm_c(PIE,settings);
% [prog, K_control, gamma, P, Z, W] = PIETOOLS_H2_control(PIE,settings);
% [prog, L_estimator, gamma, P, Z, W] = PIETOOLS_H2_estimator(PIE,settings);

% % Construct the closed-loop PIE representation, if applicable
if exist('K_control','var')
    PIE_CL = closedLoopPIE(PIE,K_control);
elseif exist('L_estimator','var')
    PIE_CL = closedLoopPIE(PIE,L_estimator,'observer');
end

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4: Simulation (See User Manual, Chapter 6 or xPIESIM/solver_PIESIM.m for more examples)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.tf=10;
opts.dt = 0.05;
opts.plot='yes';
% Declare initial conditions for PIE state
if PIE.dim<=1
    syms sx
    uinput.ic.ODE = ones(1,size(PIE.T.P,1));
    uinput.ic.PDE = ones(1,size(PIE.T.R.R0,1))*sin(sym(pi)*sx);
else
    syms sx sy
    uinput.ic.ODE = ones(1,size(PIE.T.R00,1));
    uinput.ic.PDE = ones(1,size(PIE.T.Rxx{1},1))*sin(sym(pi)*sx);
    uinput.ic.PDE = [uinput.ic.PDE,ones(1,size(PIE.T.Ryy{1},1))*sin(sym(pi)*sy)];
    uinput.ic.PDE = [uinput.ic.PDE,ones(1,size(PIE.T.R22{1},1))*sin(sym(pi)*sx)*sin(sym(pi)*sy)];
end

if ~exist('PIE_CL','var')
    % Simulate the open-loop PDE solution in PDE representation
    [solution, grid] = PIESIM(PDE,opts,uinput);
else
    % Get the number of PDE states differentiable up to each order
    ndiff = zeros(1,max(PIE.x_tab(:,end))+1);
    for jj=0:max(PIE.x_tab(:,end))
        ndiff(jj+1) = sum(PIE.x_tab(PIE.x_tab(:,end)==jj & PIE.x_tab(:,3),2));
    end
    
    % Simulate open-loop PDE solution in PIE representation
    [solution, grid] = PIESIM(PIE,opts,uinput,ndiff);
end


% Note: you can also specify time stepping and other options (opts) and user inputs, such as initial
% conditions and non-zero boundary inputs (uinput), via
% solution = PIESIM(PDE, opts),
% solution = PIESIM(PDE, uinput), or
% solution = PIESIM(PDE, opts, uinput), please consult the user's manual

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4b: Plotting solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(opts.plot,'yes') && any(ismember(fieldnames(solution),'timedep'))
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
    title('Time evolution of open-loop PDE states, $x$, plotted against space, s','Interpreter','latex');
end

% % Also simulate closed-loop response, if applicable
if exist('PIE_CL','var')
    % First, clarify that the current figures represent the open-loop
    % response
    if strcmp(opts.plot,'yes')
        figHandles_OL = findobj('Type', 'figure');
        nfigs_OL = numel(figHandles_OL);
        for jj=1:nfigs_OL
            fig_jj = figHandles_OL(jj);
            if isgraphics(fig_jj.Children(end),'Axes')
                fig_jj.Children(end).Title.String = ['Open-Loop ', fig_jj.Children(end).Title.String];
            else        %elseif isgraphics(fig_jj.Children(end),'Text')
                fig_jj.Children(end).String = ['Open-Loop ', fig_jj.Children(end).String];
            end
        end
    end

    % Update the number of states differentiable up to each order, if
    % applicable
    if exist('L_estimator','var')
        % Closed-loop system includes state estimate as well as state
        ndiff = 2*ndiff;        
        % Set initial estimate of state to 0.
        uinput.ic.ODE = [uinput.ic.ODE,0*uinput.ic.ODE];
        uinput.ic.PDE = [uinput.ic.PDE,0*uinput.ic.PDE];
    end

    % Simulate the closed-loop response
    opts.plot = 'yes';
    [solution_CL, ~] = PIESIM(PIE_CL,opts,uinput,ndiff);

    % Clarify that the new plots represent the closed-loop response
    figHandles_all = findobj('Type', 'figure');
    for jj=1:numel(figHandles_all)-nfigs_OL
        fig_jj = figHandles_all(jj);
        if exist('L_estimator','var')
            if isgraphics(fig_jj.Children(end),'Axes')
                fig_jj.Children(end).Title.String = ['True (left) and Estimated (right) ', fig_jj.Children(end).Title.String];
            else    %elseif isgraphics(fig_jj.Children(end),'Text')
                fig_jj.Children(end).String = ['True (left) and Estimated (right) ', fig_jj.Children(end).String];
            end
        else
            if isgraphics(fig_jj.Children(end),'Axes')
                fig_jj.Children(end).Title.String = ['Closed-Loop ', fig_jj.Children(end).Title.String];
            else    %elseif isgraphics(fig_jj.Children(end),'Text')
                fig_jj.Children(end).String = ['Closed-Loop ', fig_jj.Children(end).String];
            end
        end
    end

end