function PDE_out = state2pde_struct(state_in,opts)
% PDE_OUT = STATE2PDE_STRUCT(STATE_IN) takes an object of type 'state'
% and converts it to an object of type 'pde_struct', representing the same 
% terms.
%
% INPUT
% - state_in:   'state' class object, representing a vector of state
%               components, inputs, and outputs. Each component may be
%               differentiated with respect to a temporal or spatial
%               variable, and may be evaluated at a boundary of the domain.
% - opts:       Optional argument. Must be set to 'no_terms' to specify
%               that the output structure should not actually return the
%               terms in "state_in", but rather only the actual state
%               components defined in this input (i.e. states, inputs, and
%               outputs).
%
% OUTPUT
% - PDE_out:    'pde_struct' object representing the same terms as the
%               input "state_in".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 06/30/2024


% % % Process the inputs
% Check that the input state is indeed of type 'state'.
if ~isa(state_in,'state')
    error("Input must be an object of type 'state'.")
end

% Check if any options have been specified
if nargin==1
    opts = [];
else
    if ~strcmp(opts,'no_terms')
        error("Optional argument can only be set to 'no_terms'.")
    end
end

PDE_out = pde_struct();


% First, set the actual PDE variables (states, inputs, outputs), without
% declaring any terms yet.
n_states = numel(state_in);

state_type = state_in.type;        % n_states x 1 cellstr with object types
state_ID = state_in.statename;     % n_states x 1 array with object IDs
state_size = state_in.veclength;   % n_states x 1 array with object sizes
state_vars = state_in.var;         % n_states x 1 cell with variables
state_diffs = state_in.diff_order; % n_states x 1 cell with orders of derivatives


% % % Set the PDE variables in the PDE structure.

% Make sure we are not repeating objects.
[unique_IDs,idcs1,idcs2] = unique(state_ID,'stable');    % unique_IDs = IDs(idcs1);     IDs = unique_IDs(idcs2);
n_states_unique = length(unique_IDs);

% Add each unique object to the PDE structure.
comp_idcs = [0,0,0];
state_obj = cell(n_states_unique,1);       % Keep track of type of object of each state component
state_idcs = zeros(n_states_unique,1);     % Assign an index in the PDE structure to each state component
for ii=1:n_states_unique
    % Get the state number in the original list
    state_num = idcs1(ii);
    % Determine whether the state component is an actual state component,
    % an input, or an output.
    if strcmpi(state_type{state_num},'ode') || strcmpi(state_type{state_num},'pde')
        obj_ii = 'x';       obj_idx = 1;
    elseif strcmpi(state_type{state_num},'in')
        obj_ii = 'w';       obj_idx = 2;
    elseif strcmpi(state_type{state_num},'out')
        obj_ii = 'z';       obj_idx = 3;
    else
        error("Unknown state type...")
    end
    state_obj{ii} = obj_ii;

    % Declare an index for the state component.
    idx_ii = comp_idcs(obj_idx) + 1;
    comp_idcs(obj_idx) = idx_ii;
    state_idcs(ii) = idx_ii;

    % Add the state component to the PDE structure.
    PDE_out.(obj_ii){idx_ii}.ID = state_ID(state_num);
    PDE_out.(obj_ii){idx_ii}.size = state_size(state_num);
    % Set the variables of the state component.
    % For now, assume that 2 variables means one temporal variable, and one
    % spatial variable s, on domain [0,1];
    if size(state_vars{state_num},2)==1
        PDE_out.(obj_ii){idx_ii}.vars = polynomial(zeros(0,1));
        PDE_out.(obj_ii){idx_ii}.dom = zeros(0,2);
    elseif size(state_vars{state_num},2)==2
        vars_ii = state_vars{state_num};
        if ispvar(vars_ii(2))
            PDE_out.(obj_ii){idx_ii}.vars = vars_ii(2);
        else
            PDE_out.(obj_ii){idx_ii}.vars = polynomial({'s'});
        end
        PDE_out.(obj_ii){idx_ii}.dom = [0,1];
    else
        error("Conversion of 'state' objects to 'pde_struct' is currently not supported for higher-dimensional states")
    end
    % % --> take care that some variables may be evaluated at some position.
    % vars_ii = state_vars{state_num};
    % true_vars = polynomial(zeros(0,1));
    % for kk=2:size(vars_ii,2)        % start at 2, since first var is t!
    %     if ispvar(vars_ii(kk))
    %         true_vars = [true_vars; vars_ii(kk)];
    %     end
    % end
    % PDE_out.(obj_ii){idx_ii}.vars = true_vars;
    % % What about domain???
end
% % Build a table for each object to collect the ID and size.
obj_list = {'x';'y';'z';'u';'w'};
for kk=1:length(obj_list)
    obj_kk = obj_list{kk};
    obj_tab = zeros(numel(PDE_out.(obj_kk)),2);
    for ii=1:size(obj_tab,1)
        obj_tab(ii,1) = PDE_out.(obj_kk){ii}.ID;
        obj_tab(ii,2) = PDE_out.(obj_kk){ii}.size;
    end
    PDE_out.([obj_kk,'_tab']) = obj_tab;
end


if strcmp(opts,'no_terms')
    return
end

% Extend unique object types and indices to full list;
state_obj = state_obj(idcs2);
state_idcs = state_idcs(idcs2);

% % % Now, declare the actual terms.
PDE_out.free = cell(n_states,1);
for ii=1:n_states
    % Declare a single term, representing the desired derivative of the
    % state, at the desired spatial position.
    obj_ii = state_obj{ii};
    idx_ii = state_idcs(ii);
    diff_ii = state_diffs{ii};
    vars_ii = state_vars{ii};
    sz_ii = state_size(ii);

    PDE_out.free{ii}.term{1}.(obj_ii) = idx_ii;
    % Set the coefficients: just an identity matrix
    PDE_out.free{ii}.term{1}.C = eye(sz_ii);
    % Set a derivative, if applicable.
    if diff_ii(1)>0 && any(diff_ii(2:end)>0)
        error("Simultaneous differentiation with respect to temporal and spatial variable is not supported.")
    elseif diff_ii(1)>0
        % Set the desired temporal derivative of the state.
        PDE_out.free{ii}.term{1}.tdiff = diff_ii(1);
    else
        % Set the desired spatial derivative.
        PDE_out.free{ii}.term{1}.D = diff_ii(1,2:end);
    end
    % Set a location at which to evaluate, if applicable.
    PDE_out.free{ii}.term{1}.loc = vars_ii(2:end);
    % Set an empty integral.
    PDE_out.free{ii}.term{1}.I = cell(size(vars_ii,2)-1,1);

    % Keep track of the size of the free terms.
    PDE_out.free{ii}.size = sz_ii;
    % Keep track of the variables on which the free terms depend.
    true_vars = polynomial(zeros(0,1));
    for kk=2:size(vars_ii,2)
        if ispvar(vars_ii(kk))
            true_vars = [true_vars; vars_ii(kk)];
        end
    end
    PDE_out.free{ii}.vars = true_vars;
end

end