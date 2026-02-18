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
PDE = examples_PDE_library_PIETOOLS(8);

%% --- Manual Declaration Option --- 
% % To use this example, comment line 28, and
% % uncomment line 65
% pvar s t
% x = pde_var(s,[0,1]);  
% u = pde_var('control');
% z = pde_var('output',2);
% y = pde_var('sense');
% w = pde_var('input');
% eq_dyn = diff(x,t)==0*x+diff(x,s,2)+u+w;
% eq_out = z==[int(x,s,0,1);u];
% eq_obv = y==int(x,s,0,1);
% eq_bc = [subs(x,s,0)==0; subs(x,s,1)==0];
% PDE = initialize([eq_dyn;eq_out;eq_bc;eq_obv]);
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
settings.eppos = 1e-2;                % Positivity of Lyapunov Function with respect to real-valued states
settings.eppos2 = 1*1e-6;             % Positivity of Lyapunov Function with respect to spatially distributed states
settings.epneg = 0;                   % Negativity of Derivative of Lyapunov Function in both ODE and PDE state -  >0 if exponential stability desired

% % OPTIONAL: uncomment to declare SDP solver to use (defaults to solver on path)
settings.sos_opts.solver='sedumi';  % one of 'sedumi', 'mosek', 'sdpnalplus', or 'sdpt3'

% % --- Prompt for settings and choose executive automatically based on the example ---
% PIETOOLS_auto_execute

% % --- Manually run desired executives ---
% [prog,P] = lpiscript(PIE,'hinf-observer',settings);
% [prog, P] = PIETOOLS_PIE2PDEstability(PIE,settings);
% [prog, P] = PIETOOLS_PIE2PDEstability_dual(PIE,settings);
% [prog, P, gamma] = PIETOOLS_Hinf_gain(PIE,settings);
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
    uinput.ic = [ones(1,size(PIE.T.P,1)),ones(1,size(PIE.T.R.R0,1))*sin(sym(pi)*sx)];
else
    syms sx sy
    uinput.ic = ones(1,size(PIE.T.R00,1));
    uinput.ic = [uinput.ic,ones(1,size(PIE.T.Rxx{1},1))*sin(sym(pi)*sx)];
    uinput.ic = [uinput.ic,ones(1,size(PIE.T.Ryy{1},1))*sin(sym(pi)*sy)];
    uinput.ic = [uinput.ic,ones(1,size(PIE.T.R22{1},1))*sin(sym(pi)*sx)*sin(sym(pi)*sy)];
end

if ~exist('PIE_CL','var')
    % Simulate the open-loop PDE solution in PDE representation
    [solution, grid] = PIESIM(PDE,opts,uinput);
else
    % Simulate open-loop PDE solution in PIE representation
    [solution, grid] = PIESIM(PIE,opts,uinput);
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
    X = repmat(grid(:,1),1,length(t));
    n1=0;
    if PIE.dim ~=2
    pde_sol = solution.timedep.primary{2};
    else
        solution.timedep.pde{1,2,3}=solution.timedep.primary{2,3,4};
        if ~isempty(solution.timedep.pde{1})
    n=size(solution.timedep.pde{1});
    n1=n(2);
    pde_sol = solution.timedep.pde{1};
        end
        if (length(solution.timedep.pde)==2)
        if ~isempty(solution.timedep.pde{2})
    k=floor(length(grid)/2)+1;
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
        surf(t,grid(:,1),Z,'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
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
        % Set initial estimate of state to 0.
        uinput.ic.ODE = [uinput.ic.ODE,0*uinput.ic.ODE];
        uinput.ic.PDE = [uinput.ic.PDE,0*uinput.ic.PDE];
    end

    % Simulate the closed-loop response
    opts.plot = 'yes';
    [solution_CL, ~] = PIESIM(PIE_CL,opts,uinput);

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