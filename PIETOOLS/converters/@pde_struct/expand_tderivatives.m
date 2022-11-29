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
% No boundary conditions are enforced upon the newly added state
% components. Even though the new state components are temporal derivatives
% of already present components, for which BCs are given, we do not know
% how these BCs on the original states translate into BCs on the new state.
% Therefore, we will not explicitly enforce any additional BCs, instead
% letting these be implicitly enforced through the relation 
%   (d/dt) x_j = x_n.
% Nevertheless, this may introduce conservatism in e.g. stability tests.
% An updated version adding BCs may be incorporated in a later version of 
% PIETOOLS.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
nvars = size(PDE.vars,1);   % Number of spatial variables.
nx = numel(PDE.x);          % Original number of state components.
nx_new = nx;                % Updated number of state components.

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
    
    % If a higher order temporal derivative is specified, split the
    % state component as
    % \dot{x}_{ii}           = x_{n_eq+1}
    % \dot{x}_{n_eq+1}       = x_{n_eq+2}
    %      :
    % \dot{x}_{n_eq+tdiff-1} = PDE of x_{ii}
    tdiff = PDE.x{ii}.tdiff;
    vars_ii = PDE.x{ii}.vars(:,1)';
    nvars_ii = length(vars_ii);
    PDE.x = [PDE.x(:); cell(tdiff-2,1); PDE.x(ii)];
    % The added state components will have the same size and
    % variables, but will not be explicitly required to be
    % differentiable in space (this would require adding BCs).
    PDE.x_tab = [PDE.x_tab; [repmat(PDE.x_tab(ii,1:2+nvars),[tdiff-1,1]), zeros(tdiff-1,nvars)]];
    PDE.x_tab(nx_new+1:end,1) = (nx_new+1 : nx_new+tdiff-1)';
    
    % Set the equation associated to component ii to
    % \dot{x}_{ii} = x_{n_eq+1}.
    PDE.x{ii}.tdiff = 1;
    PDE.x{ii}.term = {};
    PDE.x{ii}.term{1}.x = nx_new + 1;
    PDE.x{ii}.term{1}.D = zeros(1,length(vars_ii));
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
        PDE.x{eq_num}.term{1}.D = zeros(1,length(vars_ii));
        PDE.x{eq_num}.term{1}.loc = vars_ii;
        PDE.x{eq_num}.term{1}.C = eye(PDE.x{ii}.size);
        PDE.x{eq_num}.term{1}.I = cell(nvars_ii,1);
        PDE.x{ii}.term{1}.delay = 0;
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

    % Since the current equation \dot{x}_{ii} = x_{n_eq+1} is
    % already properly specified, we just need to clean up a bit,
    % and we can move on to the next equation.
    if isfield(PDE.x{ii},'eq_num')
        PDE.x{ii} = rmfield(PDE.x{ii},'eq_num');
    end
    PDE.x{ii} = orderfields(PDE.x{ii},{'size','vars','dom','diff','tdiff','term'});

    % Update the number of state components.
    nx_new = nx_new + tdiff-1;    
    
end

% Finally, display a summary, unless otherwise indicated, and return.
PDE.has_hotd = false;
% All fields of the returned PDE should still be appropriately specified.
PDE.is_initialized = true;

if ~suppress_summary
    print_tdiff_expansion_summary(PDE,tdiff_state_tab)
end
if nx_new>nx && any(any(PDE.x_tab(nx+1:end,3:2+nvars)))
    fprintf(2,['\n Warning: No BCs have been imposed on the newly added state components representing the higher order temporal derivatives. \n',...
                 '          Results of e.g. stability tests may be very conservative.\n'])
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
    if numel(PDE.x)==1
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
    RHS = [' := ',xdot,Rcomp_idx,'(',varnames_ii_t,')'];
    
    % % % Finally, display
    MT_space = max(LHS_length_max-LHS_length,1);
    fprintf(['  ',LHS_name,repmat(' ',[1,MT_space]),RHS,'\n']);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %