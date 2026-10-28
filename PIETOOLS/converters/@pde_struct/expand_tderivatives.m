function [PDE,tdiff_state_tab] = expand_tderivatives(PDE,suppress_summary)
% PDE = expand_delays(PDE) takes a pde_struct and retruns a PDE with no
% time delays, using new state variables to represent the delayed terms.
%
% INPUTS:
% PDE:  A pde_struct class object, defining a PDE in a manner as
%       outlined in "@pde_struct/initialize".
% suppress_summary: logical value indicating whether to suppress the 
%                   output summary. Defaults to false.
% 
% OUTPUTS:
% PDE:  A pde_struct object describing the same system as the input
%       PDE, but with no higher order temporal derivative. For any PDE
%       (d^k/dt^k)x_j = F_j(x), a new state component x_n is added,
%       replacing the PDE with two new ones:
%           (d/dt) x_j = x_n
%           (d^{k-1}/dt^{k-1})x_n = F_j(x)
%       This is repeated until all PDEs are expressed in terms of first
%       order temporal derivatives.
% tdiff_state_tab:  A nx_new x 1 array of integer values, indicating for
%                   each of the new state components j=1:nx_new to which of
%                   the (original) state components it corresponds:
%                       x_j = d/dt x_{tdiff_tab(j)}.
%
% NOTES:
% For each added PDE state x_new = (d/dt)^k x_j, boundary conditions are
% imposed matching those enforced on x_j. For example, if x_j(s=0)=0, then
% also d/dt x_j(s=0)=0. Inputs in these boundary conditions are currently
% not supported, as e.g. x_j(s=0)=w would imply d/dt x_j(s=0)=d/dt w,
% requiring introduction of an input d/dt w. Boundary conditions involving
% multiple states are supported ONLY IF we have evolution equations for the
% temporal derivatives of these states as well. For example, if
% x_j(s=0)=x_i(s=0), then also d/dt x_j(s=0) = d/dt x_i(s=0), but we can
% only enforce this if we know the dynamics of d/dt x_i, i.e. if x_i is
% also constrained by a higher-order in time PDE).
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2025  PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 10/14/2022
% DJ, 01/05/2025: Account for field "ID" in object specifications.
% DJ, 06/13/2025: Update to impose BCs on newly introduced state variables;
% DJ, 10/28/2026: Skip states with tdiff=nan (no equation);

% % % Check the inputs, and initialize.
if nargin==1
    suppress_summary = false;
end

% Initialize the PDe structure, and return if there are no higher-order
% temporal derivatives.
if ~PDE.is_initialized
    PDE = initialize(PDE,true);
end
if ~PDE.has_hotd
    fprintf(['\n','No higher-order temporal derivatives were encountered.'])
    return
end

% Extract parameters.
nx = numel(PDE.x);          % Original number of state components.
nx_new = nx;                % Updated number of state components.
BC_idcs = 1:numel(PDE.BC);  % List of BCs to check.
nBCs = numel(PDE.BC);       % Original number of BCs.

% For each new state component, keep track of which of the (original) state
% components it is a temporal derivative.
tdiff_state_tab = zeros(0,1);


% % % Loop over all state components, expanding any higher-order temporal
% % % derivative (d/dt)(d^k/dt^k) x_ii = F_ii by introducing a new state 
% % % component x_jj = (d^k/dt^k) x_ii, so that d/dt x_jj = F_ii.
for ii=1:nx
    % If only a first order temporal derivative is taken, we move on.
    if ~isfield(PDE.x{ii},'tdiff') || PDE.x{ii}.tdiff==1
        continue
    end
    % A state without equation has no temporal derivative to expand.
    if isnan(PDE.x{ii}.tdiff)                                               % DJ, 10/28/2026
        continue                                                            
    end 
    
    % If a higher order temporal derivative is specified, split the
    % state component as
    % \dot{x}_{ii}           = x_{n_eq+1}
    % \dot{x}_{n_eq+1}       = x_{n_eq+2}
    %      :
    % \dot{x}_{n_eq+tdiff-1} = PDE of x_{ii}
    tdiff = PDE.x{ii}.tdiff;
    vars_ii = PDE.x{ii}.vars(:,1)';
    nvars_ii = numel(vars_ii);
    PDE.x = [PDE.x(:); cell(tdiff-2,1); PDE.x(ii)];
    % The added state components will have the same size and
    % variables, and be assumed to be differentiable up to the same order.
    PDE.x_tab = [PDE.x_tab; [repmat(PDE.x_tab(ii,1:end),[tdiff-1,1])]];
    PDE.x_tab(nx_new+1:end,1) = (nx_new+1 : nx_new+tdiff-1)';
    
    % Set the equation associated to component ii to
    % \dot{x}_{ii} = x_{n_eq+1}.
    PDE.x{ii}.tdiff = 1;
    PDE.x{ii}.term = {};
    PDE.x{ii}.term{1}.x = nx_new + 1;
    PDE.x{ii}.term{1}.D = zeros(1,nvars_ii);
    PDE.x{ii}.term{1}.loc = vars_ii;
    PDE.x{ii}.term{1}.C = eye(PDE.x{ii}.size);
    PDE.x{ii}.term{1}.I = cell(nvars_ii,1);
    PDE.x{ii}.term{1}.delay = 0;
    
    % Set the equations associated to the other components as
    % \dot{x}_{n_eq+kk} = x_{n_eq+kk+1}.
    for extra_eq = 1:tdiff-2
        eq_num = nx_new + extra_eq;
        PDE.x{eq_num}.eq_num = ii;
        PDE.x{eq_num}.vars = PDE.x{ii}.vars;
        PDE.x{eq_num}.dom = PDE.x{ii}.dom;
        PDE.x{eq_num}.diff = zeros(size(PDE.x{ii}.diff));
        PDE.x{eq_num}.tdiff = 1;
        
        PDE.x{eq_num}.term{1}.x = eq_num + 1;
        PDE.x{eq_num}.term{1}.D = zeros(1,nvars_ii);
        PDE.x{eq_num}.term{1}.loc = vars_ii;
        PDE.x{eq_num}.term{1}.C = eye(PDE.x{ii}.size);
        PDE.x{eq_num}.term{1}.I = cell(nvars_ii,1);
        PDE.x{eq_num}.term{1}.delay = 0;
    end
    % For the final component, we already have the equation, 
    % \dot{x}_{n_eq+tdiff-1} = F_ii(x) where (d/dt)^tdiff x_ii = F_ii(x).
    % We just need to update the order of the temporal derivative,
    % and the order of differentiability in the spatial variables.
    PDE.x{end}.tdiff = 1;
    PDE.x{end}.diff = zeros(size(PDE.x{ii}.diff));
    
    % Keep track of which state component each new state component is the
    % temporal derivative.
    tdiff_state_tab = [tdiff_state_tab; [ii; nx_new+(1:tdiff-2)']];

    % Declare boundary conditions for newly introduced state variables
    for BC_num = BC_idcs                                                    % DJ, 06/13/2025
        % Check whether any of the terms involve state variable ii
        x_idcs = zeros(1,numel(PDE.BC{BC_num}.term));
        for jj=1:numel(PDE.BC{BC_num}.term)
            if isfield(PDE.BC{BC_num}.term{jj},'x')
                x_idcs(jj) = PDE.BC{BC_num}.term{jj}.x;
            end
        end
        if ~any(x_idcs==ii)
            continue
        end
        % Check if boundary condition already correponds to temporal
        % derivative of other state
        if ~isfield(PDE.BC{BC_num},'tdiff')
            % The BC is one already specified by the user
            % --> impose the same boundary condition on the temporal
            % derivative of the state variable
            BC_idcs = setdiff(BC_idcs,BC_num);
            
            % Make sure the BC doesn't involves any other inputs
            if ~all(x_idcs)
                error("Expansion of temporal derivatives for states with boundary inputs is currently not supported.")
            end
            % Impose same boundary condition on each state diff(x,'t',kk)
            BCs_tmp = repmat(PDE.BC(BC_num),[tdiff-1,1]);
            BC_tab_tmp = [PDE.BC_tab(end,1)+(1:tdiff-1)',repmat(PDE.BC_tab(BC_num,2:end),[tdiff-1,1])];
            for kk=1:tdiff-1
                for jj=find(x_idcs==ii)
                    % Replace state index with that of diff(x_{ii},'t',kk)
                    BCs_tmp{kk}.term{jj}.x = nx_new+kk;
                end
            end
            % Keep track of any other state components the BC involves
            if any(x_idcs~=ii)
                % BC involves other state x_{jj}
                % --> BC on temporal derivative of x_{ii} will involve temporal
                %       derivative of x_{jj}
                % We will need to modify the BC to account for this
                BC_idcs = [BC_idcs,numel(PDE.BC)+(1:tdiff-1)];
                for kk=1:tdiff-1
                    BCs_tmp{kk}.tdiff = kk;
                    BCs_tmp{kk}.tdiff_idcs = find(x_idcs~=ii);
                end
            end        
            PDE.BC = [PDE.BC; BCs_tmp];
            PDE.BC_tab = [PDE.BC_tab; BC_tab_tmp];
        else
            % The BC already corresponds to a temporal derivative of an 
            % earlier state component.
            % --> modify index of current state component to match this
            % temporal derivative.

            % We can't take a higher-order temporal derivative than tdiff
            if PDE.BC{BC_num}.tdiff>tdiff-1
                error("Expansion of temporal derivatives is not supported: boundary conditions on temporal derivative of PDE state are ambiguous.")
            end
            trm_idcs = PDE.BC{BC_num}.tdiff_idcs;
            for jj=trm_idcs(x_idcs(trm_idcs)==ii)
                PDE.BC{BC_num}.term{jj}.x = nx_new + PDE.BC{BC_num}.tdiff;
            end
            % Remove state ii from list of indices that need to be updated.
            trm_idcs = trm_idcs(x_idcs(trm_idcs)~=ii);
            if isempty(trm_idcs)
                % All state indices in the BC have been updated.
                BC_idcs = setdiff(BC_idcs,BC_num);
                PDE.BC{BC_num} = rmfield(PDE.BC{BC_num},{'tdiff','tdiff_idcs'});
            else
                % There are still state indices that need to be updated.
                PDE.BC{BC_num}.tdiff_idcs = trm_idcs;
            end
        end
    end

    % Since the current equation \dot{x}_{ii} = x_{n_eq+1} is
    % already properly specified, we just need to clean up a bit,
    % and we can move on to the next equation.
    if isfield(PDE.x{ii},'eq_num')
        PDE.x{ii} = rmfield(PDE.x{ii},'eq_num');
    end
    PDE.x{ii} = orderfields(PDE.x{ii},{'ID','size','vars','dom','diff','tdiff','term'});        % DJ, 01/05/2025

    % Update the number of state components.
    nx_new = nx_new + tdiff-1;    
    
end

% Check if there are any "unfinished" BCs (e.g. x_j(s=0) = x_i(s=0) implies
% d/dt x_j(s=0) = d/dt x_i(s=0) but we have no state variable representing 
% d/dt x_i)
if any(BC_idcs>nBCs)
    warning("Boundary conditions on temporal derivative of certain PDE state are ambiguous; returned system may not be well-posed.")
end
% Get rid of BCs on auxiliary variables which cannot be completed.
PDE.BC(BC_idcs(BC_idcs>nBCs)) = [];
PDE.BC_tab(BC_idcs(BC_idcs>nBCs),:) = [];

% Finally, display a summary, unless otherwise indicated, and return.
PDE.has_hotd = false;
% All fields of the returned PDE should still be appropriately specified.
PDE.is_initialized = true;

if ~suppress_summary
    print_tdiff_expansion_summary(PDE,tdiff_state_tab)
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function print_tdiff_expansion_summary(PDE,tdiff_state_tab)
% print_tdifff_expanion_summary(PDE,obj,ncomps_x)
% prints in the command window some information concerning how many new
% state components have been introduced to represent the higher-order
% temporal derivatives, and which of the original state components these
% correspond to.
%
% INPUTS:
% - PDE:    A "pde_struct" class object defining a PDE.
% - tdiff_state_tab:    A nx_new x 1 array of integers indicating for each
%                       of the j=1:nx_new added state components of which
%                       state component they are the temporal derivative:
%                       x_j = d/dt x_{tdiff_tab(j)}.
%
% OUTPUTS:
% Displays information in the command window concerning what state
% components have been added to represent the higher order temporal
% derivatives, and to what state components they correspond.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nx_full = numel(PDE.x);
nx_new = numel(tdiff_state_tab);
nx = nx_full - nx_new;
if nx_new==0
    fprintf(['\n','No higher-order temporal derivatives were encountered.'])
elseif nx_new==1
    fprintf(['\n','Added 1 state component: \n']);
else
    fprintf(['\n','Added ',num2str(nx_new),' state components: \n']);
end


% Use UNICODE to add subscript indices to different components.
sub_num = mat2cell([repmat('\x208',[10,1]),num2str((0:9)')],ones(10,1),6);
partial = '\x2202';
sub_t = '\x209C';
%int = '\x222B';
%xdot = '\x1E8B';
xdot = [partial,sub_t,' x'];
% Determine how many subscripts will be needed.
n_digits = 1+floor(log(nx_full)/log(10));


% % For the purpose of clearer display, we make an estimate of how long the
% % display will be. 
% First, determine which variables appear in the different state variables.
nvars = size(PDE.vars,1);
global_vars_obj = PDE.vars(any(PDE.x_tab(:,3:2+nvars),1),:);
if isempty(global_vars_obj)
    global_varnames = '(t)';
else
    global_varnames =  global_vars_obj(1,1).varname{1};
    for kk=2:size(PDE.vars,1)
        global_varnames = [global_varnames,',',PDE.vars(kk,1).varname{1}];
    end
    global_varnames = ['(t,',global_varnames,')'];
end
% Then, estimate the size of the string of characters denoting the components
% "PDE.x{ii}".
nvars_max = max(sum(PDE.x_tab(:,3:2+nvars),2));
lngth_varnames_mean = ceil(length(global_varnames)*nvars_max/nvars);
LHS_length_max = 1+1 + n_digits + lngth_varnames_mean+2; % e.g. ' x13(t,s1,s2,s3)', 


% % For each of the components, display its size, and which variables it
% % depends on.
for ii=1:nx_new
    Lstate_num = ii + nx;   % State index associated to the new state.
    comp_ii = PDE.x{Lstate_num};
    
    % Establish the names of the variables on which the component depends.
    if isempty(comp_ii.vars)
        varnames_ii = '';
        varnames_ii_t = 't';    % All components may vary in time
    else
        % Set a comma separated list.
        varnames_ii = comp_ii.vars(1,1).varname{1};
        for kk=2:size(comp_ii.vars,1)
            varnames_ii = [varnames_ii,',',comp_ii.vars(kk,1).varname{1}];
        end
        varnames_ii_t = ['t,',varnames_ii]; % All components may vary in time
    end
    
    % Establish the (subscript) index for the component.
    LHS_length = length(varnames_ii_t);
    if isscalar(PDE.x)
        % There is only one state component --> no need to give index
        Lcomp_idx = '';
    elseif numel(PDE.x)<=9
        % The state number consists of a single decimal.
        Lcomp_idx = sub_num{Lstate_num+1};
        LHS_length = LHS_length + 1;
    else
        % The state number consists of multiple decimals.
        Lcomp_idx = cell2mat(sub_num(str2num(num2str(Lstate_num)')+1)');
        LHS_length = LHS_length + length(num2str(Lstate_num));
    end
    % Set the name of the component, including its depdence on spatial
    % variables.
    LHS_name = [' x',Lcomp_idx,'(',varnames_ii_t,')'];
    LHS_length = 1 + LHS_length + 3;
        
    % For the added state components, indicate of which state component
    % they are the temporal derivative.
    Rstate_num = tdiff_state_tab(ii);
    Rcomp_idx = cell2mat(sub_num(str2num(num2str(Rstate_num)')+1)');
    RHS = [' :=  ',xdot,Rcomp_idx,'(',varnames_ii_t,')'];
    
    % % % Finally, display
    MT_space = max(LHS_length_max-LHS_length,1);
    fprintf(['  ',LHS_name,repmat(' ',[1,MT_space]),RHS,'\n']);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %