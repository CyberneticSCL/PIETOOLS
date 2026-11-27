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

% This script declares a linear coupled ODE-PDE (Step 1), converts it to
% a PIE (Step 2), runs the executive uncommented in Step 3 (none by
% default), and simulates it with PIESIM (Step 4).
%
% INPUT
% - none; Step 1 takes the PDE from the example library (example 17), the
%   GUI, or the manual declaration, whichever is uncommented.
%
% OUTPUT
% - none; the workspace receives the PDE, its PIE, the LPI settings, the
%   outputs of the executive uncommented in Step 3 (none by default; with
%   PIE_CL, the closed-loop PIE, for a controller or estimator), and the
%   PIESIM solution of Step 4.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2026 PIETOOLS Team
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
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding PIETOOLS Team (PIETOOLS 2021b; in the repository since 12/31/2021)
% YP, MMP, DB, 11/21/2026: Changes since PIETOOLS 2025: Step 1 loads library
%   example 17; Step 3 takes the 'heavy' settings of the container
%   executives and lists these executives; Step 4 passes the initial
%   conditions as one vector uinput.ic, sized from PIE.space (for PIESIM of
%   a pde_struct, in the order of PDE.x_tab), with a zero initial estimate
%   per kind of space for an estimator, and calls PIESIM without ndiff;
%   Step 4b plots the PIESIM output fields primary and grid, in 2-D the
%   states in s1 and in (s1,s2) at the middle node in s2; help and comments
%   revised.
% MMP, 11/27/2026: copvar, cdopvar renamed opvar, dopvar, with the
%   functions of the rename table; names only.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Declare a PDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % --- GUI option ---
% PIETOOLS_PDE_GUI
%%
% % --- Example Library Option ---
% The example library contains 72 examples of ODE-PDE systems,
% each with an associated default executive to run. Call the library with
% an index 1 through 72 as argument to extract that example, and run the
% associated executive, if desired.
PDE = examples_PDE_library_PIETOOLS(17);

%% --- Manual Declaration Option --- 
% % To use this example, comment out the library call above, and
% % uncomment the lines below
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
% settings = lpisettings('heavy');                                          % YP, MMP, DB, 11/21/2026 (was)
settings = lpisettings_sop('heavy');  % tiers of the container executives   % YP, MMP, DB, 11/21/2026
settings.eppos = 1e-2;                % Positivity of Lyapunov Function with respect to real-valued states
settings.eppos2 = 1*1e-6;             % Positivity of Lyapunov Function with respect to spatially distributed states
settings.epneg = 0;                   % Negativity of Derivative of Lyapunov Function in both ODE and PDE state -  >0 if exponential stability desired

% % OPTIONAL: uncomment to declare SDP solver to use (defaults to solver on path)
settings.sos_opts.solver='sedumi';  % one of 'sedumi', 'mosek', 'sdpnalplus', or 'sdpt3'

% % --- Prompt for settings and choose executive automatically based on the example ---
% PIETOOLS_auto_execute                                                     % YP, MMP, DB, 11/21/2026 (was)
% PIETOOLS_auto_execute_sop                                                 % YP, MMP, DB, 11/21/2026

% % --- Manually run desired executives ---
% The container executives (PIETOOLS_*_sop) take the opvar PIE that        % YP, MMP, DB, 11/21/2026
% convert_PIETOOLS_PDE returns and return the outputs named, then an info   % YP, MMP, DB, 11/21/2026
% struct; the third argument is the form ('Q', 'Qd': PIE2PDE forms; 'o',    % YP, MMP, DB, 11/21/2026
% 'c': observability and controllability H2 norms).                         % YP, MMP, DB, 11/21/2026
% [prog,P] = lpiscript(PIE,'hinf-observer',settings);                       % YP, MMP, DB, 11/21/2026 (was)
% [prog,P] = lpiscript_sop(PIE,'hinf-observer',settings);                   % YP, MMP, DB, 11/21/2026
% [prog, P] = PIETOOLS_PIE2PDEstability(PIE,settings);                      % YP, MMP, DB, 11/21/2026 (was)
% [prog, P] = PIETOOLS_stability_sop(PIE,settings,'Q');                     % YP, MMP, DB, 11/21/2026
% [prog, P] = PIETOOLS_PIE2PDEstability_dual(PIE,settings);                 % YP, MMP, DB, 11/21/2026 (was)
% [prog, P] = PIETOOLS_stability_sop(PIE,settings,'Qd');                    % YP, MMP, DB, 11/21/2026
% [prog, P, gamma] = PIETOOLS_Hinf_gain(PIE,settings);                      % YP, MMP, DB, 11/21/2026 (was)
% [prog, P, gamma] = PIETOOLS_Hinf_gain_sop(PIE,settings,'Q');              % YP, MMP, DB, 11/21/2026
% [prog, P, gamma] = PIETOOLS_Hinf_gain_dual(PIE,settings);                 % YP, MMP, DB, 11/21/2026 (was)
% [prog, P, gamma] = PIETOOLS_Hinf_gain_sop(PIE,settings,'Qd');             % YP, MMP, DB, 11/21/2026
% [prog, K_control, gamma, P, Z] = PIETOOLS_Hinf_control(PIE,settings);     % YP, MMP, DB, 11/21/2026 (was)
% [prog, K_control, gamma, P, Z] = PIETOOLS_Hinf_control_sop(PIE,settings); % YP, MMP, DB, 11/21/2026
% [prog, L_estimator, gamma, P, Z] = PIETOOLS_Hinf_estimator(PIE,settings); % YP, MMP, DB, 11/21/2026 (was)
% [prog, L_estimator, gamma, P, Z] = PIETOOLS_Hinf_estimator_sop(PIE,settings); % YP, MMP, DB, 11/21/2026
% [prog, Wo, gamma] = PIETOOLS_H2_norm_o(PIE,settings);                     % YP, MMP, DB, 11/21/2026 (was)
% [prog, Wo, gamma] = PIETOOLS_H2_norm_sop(PIE,settings,'o');               % YP, MMP, DB, 11/21/2026
% [prog, Wc, gamma] = PIETOOLS_H2_norm_c(PIE,settings);                     % YP, MMP, DB, 11/21/2026 (was)
% [prog, Wc, gamma] = PIETOOLS_H2_norm_sop(PIE,settings,'c');               % YP, MMP, DB, 11/21/2026
% [prog, K_control, gamma, P, Z, W] = PIETOOLS_H2_control(PIE,settings);    % YP, MMP, DB, 11/21/2026 (was)
% [prog, K_control, gamma, P, Z, W] = PIETOOLS_H2_control_sop(PIE,settings); % YP, MMP, DB, 11/21/2026
% [prog, L_estimator, gamma, P, Z, W] = PIETOOLS_H2_estimator(PIE,settings); % YP, MMP, DB, 11/21/2026 (was)
% [prog, L_estimator, gamma, P, Z, W] = PIETOOLS_H2_estimator_sop(PIE,settings); % YP, MMP, DB, 11/21/2026

% % Construct the closed-loop PIE representation, if applicable
if exist('K_control','var')
    PIE_CL = closedLoopPIE(PIE,K_control);
elseif exist('L_estimator','var')
    PIE_CL = closedLoopPIE(PIE,L_estimator,'observer');
end

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4: Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final time, time step, and plotting by PIESIM
opts.tf=10;
opts.dt = 0.05;
opts.plot='yes';
% Declare initial conditions for PIE state
% Component sizes and spatial dependence of the PIE state, from PIE.space;  % YP, MMP, DB, 11/21/2026
% an uninitialized PIE (legacy terms format) is initialized in a copy       % YP, MMP, DB, 11/21/2026
PIE_init = PIE;                                                             % YP, MMP, DB, 11/21/2026
if ~PIE_init.is_initialized,    PIE_init = initialize(PIE_init);    end     % YP, MMP, DB, 11/21/2026
xdims = PIE_init.space.x.dims;  xdep = PIE_init.space.x.depmat;             % YP, MMP, DB, 11/21/2026
if PIE.dim<=1
    syms sx
    uinput.ic = [ones(1,sum(xdims(~any(xdep,2)))),ones(1,sum(xdims(any(xdep,2))))*sin(sym(pi)*sx)]; % YP, MMP, DB, 11/21/2026
else
    syms sx sy
    uinput.ic = ones(1,sum(xdims(~any(xdep,2))));                           % YP, MMP, DB, 11/21/2026
    uinput.ic = [uinput.ic,ones(1,sum(xdims(xdep(:,1) & ~xdep(:,2))))*sin(sym(pi)*sx)]; % YP, MMP, DB, 11/21/2026
    uinput.ic = [uinput.ic,ones(1,sum(xdims(~xdep(:,1) & xdep(:,2))))*sin(sym(pi)*sy)]; % YP, MMP, DB, 11/21/2026
    uinput.ic = [uinput.ic,ones(1,sum(xdims(xdep(:,1) & xdep(:,2))))*sin(sym(pi)*sx)*sin(sym(pi)*sy)]; % YP, MMP, DB, 11/21/2026
end
% PIESIM(PDE,...) below reads one entry per scalar state component of the   % YP, MMP, DB, 11/21/2026
% PDE, in the order of PDE.x_tab (size in column 2, dependence on each      % YP, MMP, DB, 11/21/2026
% variable in columns 3:2+dim); the other PDE formats and the PIE keep the  % YP, MMP, DB, 11/21/2026
% layout above, by kind                                                     % YP, MMP, DB, 11/21/2026
if isa(PDE,'pde_struct') && ~exist('PIE_CL','var')                          % YP, MMP, DB, 11/21/2026
    PDE_init = PDE;                                                         % YP, MMP, DB, 11/21/2026
    if ~PDE_init.is_initialized,    PDE_init = initialize(PDE_init,true);   end % YP, MMP, DB, 11/21/2026
    syms sx sy                                                              % YP, MMP, DB, 11/21/2026
    % The entry of each kind 1-4: R, L2[s1], L2[s2], L2[s1,s2]              % YP, MMP, DB, 11/21/2026
    ic_kind = [sym(1),sin(sym(pi)*sx),sin(sym(pi)*sy),sin(sym(pi)*sx)*sin(sym(pi)*sy)]; % YP, MMP, DB, 11/21/2026
    nd = min(PDE_init.dim,2);   dep = false(size(PDE_init.x_tab,1),2);      % YP, MMP, DB, 11/21/2026
    dep(:,1:nd) = PDE_init.x_tab(:,3:2+nd)~=0;                              % YP, MMP, DB, 11/21/2026
    uinput.ic = ic_kind(repelem(1+dep(:,1)+2*dep(:,2),PDE_init.x_tab(:,2))); % YP, MMP, DB, 11/21/2026
end                                                                         % YP, MMP, DB, 11/21/2026

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
% solution = PIESIM(PDE, opts, uinput).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4b: Plotting solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the distributed states against t and s when PIESIM did not plot them
if ~strcmp(opts.plot,'yes') && any(ismember(fieldnames(solution),'timedep'))
    % Time grid
    t = solution.timedep.dtime;
    X = repmat(grid(:,1),1,length(t));
    n1=0;
    % Distributed states: in 1-D primary{2}; in 2-D the states in s1, then
    % the states in (s1,s2) at the middle grid point in s2
    if PIE.dim ~=2
    pde_sol = solution.timedep.primary{2};
    else
        solution.timedep.pde=solution.timedep.primary([2,4]);  % s1, 2-D    % YP, MMP, DB, 11/21/2026
        if ~isempty(solution.timedep.pde{1})
    n=size(solution.timedep.pde{1});
    n1=n(2);
    pde_sol = solution.timedep.pde{1};
        end
        if (length(solution.timedep.pde)==2)
        if ~isempty(solution.timedep.pde{2})
    k=floor(length(grid{2})/2)+1;  % middle node of the s2 grid             % YP, MMP, DB, 11/21/2026
    n=size(solution.timedep.pde{2});
    pde_sol(1:n(1),n1+1:n1+n(3),1:n(4)) = solution.timedep.pde{2}(1:n(1),k,1:n(3),1:n(4));
        end
        end
    end
    if iscell(grid), sgrid=grid{1}; else, sgrid=grid(:,1); end  % s1 grid   % YP, MMP, DB, 11/21/2026
    % One surface plot over (t,s) per state component
    t = repmat(t,size(pde_sol,1),1);
    figure;
    for i=1:size(pde_sol,2)
        Z = squeeze(pde_sol(:,i,:));
        ax(i)=subplot(size(pde_sol,2),1,i);
        surf(t,sgrid,Z,'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row'); % YP, MMP, DB, 11/21/2026
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
%       uinput.ic.ODE = [uinput.ic.ODE,0*uinput.ic.ODE];                    % YP, MMP, DB, 11/21/2026 (was)
%       uinput.ic.PDE = [uinput.ic.PDE,0*uinput.ic.PDE];                    % YP, MMP, DB, 11/21/2026 (was)
        % uinput.ic is one vector, the kinds of space in turn (Step 4); in  % YP, MMP, DB, 11/21/2026
        % each kind the closed-loop state [x;xhat] has x, then xhat.        % YP, MMP, DB, 11/21/2026
        kinds = 1 + any(xdep,2);                            % R, L2         % YP, MMP, DB, 11/21/2026
        if PIE.dim>1,   kinds = 1 + xdep(:,1) + 2*xdep(:,2);    end         % YP, MMP, DB, 11/21/2026
        ick = mat2cell(uinput.ic,1,accumarray(kinds,xdims,[2+2*(PIE.dim>1),1])'); % YP, MMP, DB, 11/21/2026
        uinput.ic = [ick; cellfun(@(c) 0*c,ick,'UniformOutput',false)];     % YP, MMP, DB, 11/21/2026
        uinput.ic = [uinput.ic{:}];                                         % YP, MMP, DB, 11/21/2026
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
