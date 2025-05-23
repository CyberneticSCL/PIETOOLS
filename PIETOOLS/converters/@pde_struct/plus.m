function [PDE_out] = plus(PDE_1,PDE_2)
% PDE_OUT = PLUS(PDE_1,PDE_2) builds a PDE representing the sum of the 
% terms in the PDE objects PDE_1 and PDE_2.
%
% INPUT
% - PDE_1:      'pde_struct' object representing either a set of equations
%               ([d/dt x=... ; y=...; z=...; 0=...]), or a set of free 
%               terms to be used to declare an equation.
% - PDE_2:      'pde_struct' object representing a set of free terms to 
%               add to the equations or terms specified by PDE_1. 
%
% OUTPUT
% - PDE_out:    'pde_struct' object representing either the same equations
%               as in PDE_1 but now with the terms from PDE_2 added, or a
%               new set of free terms (collected in PDE_out.free)
%               corresponding to the sum of the terms in PDE_1 and PDE_2;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024  PIETOOLS Team
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
% DJ, 06/23/2024: Initial coding
% DJ, 01/03/2025: Update to assume a loose PDE variable is specified as a
%                   single free term, see also update to "pde_var".
%                   Also account for added "is_zero" field;
% DJ, 01/04/2025: Bugfix for empty terms;


% % % Process the inputs
% % Check that inputs are of appropriate type.
if isa(PDE_1,'state')
    PDE_1 = state2pde_struct(PDE_1);
elseif isa(PDE_1,'terms')
    PDE_1 = terms2pde_struct(PDE_1);
end
if ~isa(PDE_1,'pde_struct') 
    % Summation with non-pde_struct objects is supported only for zeros.
    if isa(PDE_1,'polynomial') && isdouble(PDE_1)
        PDE_1 = double(PDE_1);
    end
    if isnumeric(PDE_1) && all(all(PDE_1==0))
        % Summation with zero just returns the input terms.
        PDE_out = PDE_2;
        return
    else
        error("Summation of 'pde_struct' objects with non-'pde_struct' objects is not supported.")
    end
end
if isa(PDE_2,'state')
    PDE_2 = state2pde_struct(PDE_2);
elseif isa(PDE_2,'terms')
    PDE_2 = terms2pde_struct(PDE_2);
end
if ~isa(PDE_2,'pde_struct') 
    % Summation with non-pde_struct objects is supported only for zeros.
    if isa(PDE_2,'polynomial') && isdouble(PDE_2)
        PDE_2 = double(PDE_2);
    end
    if isnumeric(PDE_2) && all(all(PDE_2==0))
        % Summation with zero just returns the input terms.
        PDE_out = PDE_1;
        return
    else
        error("Summation of 'pde_struct' objects with non-'pde_struct' objects is not supported.")
    end
end

% % Make sure at least one of the inputs correspond to loose terms (we
% % can't add full PDEs
if ~is_pde_term(PDE_1) && ~is_pde_term(PDE_2)
    error("Addition of full PDE equations is not supported.")
elseif ~is_pde_term(PDE_2)
    % Trying to add terms on left-hand side of PDE
    % --> move to right-hand side
    PDE_out = plus(PDE_2,-1*PDE_1);
end

% % Determine what type of objects we can extract from PDE_1
if is_pde_term(PDE_1)
    objs_1 = {'free'};
    n_eqs_arr = numel(PDE_1.free);
    n_eqs_1 = n_eqs_arr;
    ncomps_1 = numel(PDE_1.free);
else
    % Assume order [d/dt x; y; z; BC]
    objs_1 = {'x';'y';'z';'BC'};
    ncomps_arr = [numel(PDE_1.x); numel(PDE_1.y); numel(PDE_1.z); numel(PDE_1.BC)];
    ncomps_1 = sum(ncomps_arr);
    % Check for how many objects a PDE has actually been specified.
    n_eqs_arr = ncomps_arr;
    include_idcs_1 = cell(1,length(ncomps_arr));
    idx = 0;
    for kk=1:numel(objs_1)
        include_idcs_1{kk} = [];
        for ii=1:numel(PDE_1.(objs_1{kk}))
            idx = idx+1;
            if ~isfield(PDE_1.(objs_1{kk}){ii},'term') || isempty(PDE_1.(objs_1{kk}){ii}.term)
                % No equation has been specified for this object.
                n_eqs_arr(kk) = n_eqs_arr(kk)-1;
            else
                include_idcs_1{kk} = [include_idcs_1{kk};ii];
            end
        end
    end
    n_eqs_1 = sum(n_eqs_arr);
end
ncomps_2 = numel(PDE_2.free);
nn_eqs_arr = cumsum([0;n_eqs_arr]);   % If nncomps_arr(j-1)<ii<=nncomps_arr(j), then equation ii corresponds to object objs_1{j};

% Assume terms need to be added only to components for which an equation
% has been specified.
add_eqs_only = true;
if n_eqs_1~=ncomps_2
    if ncomps_1~=ncomps_2
        error("The numbers of rows of terms to add do not match.")
    else
        % assume terms need to be added to all components.
        nn_eqs_arr = cumsum([0;ncomps_arr]);
        add_eqs_only = false;
    end
end

% % % Perform the actual addition
% % First express the PDEs in terms of the same components
[PDE_1_out,PDE_2_out,new_obnums_1] = pde_common_basis(PDE_1,PDE_2);

% % Initialize the output as the first PDE
PDE_out = PDE_1_out;

% % Loop over all equations, and add the terms.
if strcmp(objs_1,'free')
    % Both PDE_1 and PDE_2 correspond to loose terms, yet to define an
    % actual equation.
    for ii=1:ncomps_2
        % Check that the number of rows in the terms to add make sense.
        if isfield(PDE_1_out.free{ii},'size') && isfield(PDE_2_out.free{ii},'size') &&...
                PDE_1_out.free{ii}.size~=PDE_2_out.free{ii}.size
            error("The terms to add have different numbers of rows.")
        end
        % Keep track of which variables appear in the terms.
        if isfield(PDE_1_out.free{ii},'vars') && ~isempty(PDE_1_out.free{ii}.vars)
            varnames_1 = PDE_1_out.free{ii}.vars.varname;
        else
            varnames_1 = {};
        end
        if isfield(PDE_2_out.free{ii},'vars') && ~isempty(PDE_2_out.free{ii}.vars)
            varnames_2 = PDE_2_out.free{ii}.vars.varname;
        else
            varnames_2 = {};
        end
        PDE_out.free{ii}.vars = polynomial(unique([varnames_1;varnames_2]));

        % Add the actual terms.                                             % DJ, 01/04/2025
        if isfield(PDE_1_out.free{ii},'term') && ~isempty(PDE_1_out.free{ii}.term) &&...
                isfield(PDE_2_out.free{ii},'term') && ~isempty(PDE_2_out.free{ii}.term)
            PDE_out.free{ii}.term = [PDE_1_out.free{ii}.term, PDE_2_out.free{ii}.term];
        elseif isfield(PDE_1_out.free{ii},'term') && ~isempty(PDE_1_out.free{ii}.term)
            PDE_out.free{ii}.term = PDE_1_out.free{ii}.term;
            continue
        elseif isfield(PDE_2_out.free{ii},'term') && ~isempty(PDE_2_out.free{ii}.term)
            PDE_out.free{ii}.term = PDE_2_out.free{ii}.term;
            continue
        else
            PDE_out.free{ii}.term = {};
            continue
        end
        % Check that the resulting equation does not contain multiple
        % left-hand side objects (temporal derivatives of states, outputs).
        % These objects should always be moved to the first term.
        [is_LHS_2,obj_2,eq_num_2,tdiff_2] = is_LHS_term(PDE_2_out.free{ii}.term{1});
        if is_LHS_2
            [is_LHS_1,obj_1,eq_num_1,tdiff_1] = is_LHS_term(PDE_1_out.free{ii}.term{1});
            if is_LHS_1
                % Both PDEs contain a ``left-hand side object''. This is
                % only allowed if they cancel
                if strcmp(obj_1,obj_2) && eq_num_1==eq_num_2 && tdiff_1==tdiff_2 &&...
                       all(all(isequal(PDE_1_out.free{ii}.term{1}.C-PDE_2_out.free{ii}.term{1}.C,0)))
                    PDE_out.free{ii}.term = PDE_out.free{ii}.term(2:end);
                else
                    error("Any PDE equation can contain at most 1 temporal derivative of a state, or 1 output signal.")
                end
            else
                % PDE_2 contains a ``left-hand side object''. Move this
                % term to the start.
                PDE_out.free{ii}.term = [PDE_2_out.free{ii}.term(1),PDE_1_out.free{ii}.term,PDE_2_out.free{ii}.term(2:end)];
            end
        end
    end
else
    % PDE_1 corresponds to an actual (completed) set of PDE equations to
    % which to add the terms in PDE_2
    for ii=1:ncomps_2
        % Make sure PDE_2 does not include any ``left-hand side object'',
        % i.e. a temporal derivative of a state variable, or an output
        % variable, which could certainly not be added to an already
        % completed PDE.
        if ~isfield(PDE_2.free{ii},'term') || isempty(PDE_2.free{ii}.term)
            continue
        elseif is_LHS_term(PDE_2.free{ii}.term{1})
            error("Any PDE equation can contain at most 1 temporal derivative of a state, or 1 output signal.")
        end

        % Establish which type of object corresponds to equation ii in
        % PDE_1 (not in PDE_1_out!)
        ob_idx = find(ii<=nn_eqs_arr,1,'first')-1;
        obj_1 = objs_1{ob_idx};
        eq_num_1 = ii - nn_eqs_arr(ob_idx);
        if add_eqs_only
            % Adjust equation number in case only already specified
            % equations get new terms.
            eq_num_1 = include_idcs_1{ob_idx}(eq_num_1);
        end

        % Check that the size of the terms in PDE_1 and PDE_2 match.
        if PDE_1.(obj_1){eq_num_1}.size~=PDE_2.free{ii}.size
            error("The row-dimension of the terms to add doesn't match.")
        end
        % Check that PDE_2 does not introduce variables on which the object
        % in PDE_1 does not depend
        if any(~ismember(PDE_2.free{ii}.vars.varname,PDE_1.(obj_1){eq_num_1}.vars.varname))
            error("The proposed terms cannot be added to desired PDE equations: the spatial variables do not match.")
        end
        % Establish which equation in the output PDE corresponds to
        % equation ii in PDE_1
        if strcmp(obj_1,'x')
            eq_num = new_obnums_1{1}(eq_num_1);
        elseif strcmp(obj_1,'y')
            eq_num = new_obnums_1{2}(eq_num_1);
        elseif strcmp(obj_1,'z')
            eq_num = new_obnums_1{3}(eq_num_1);
        else
            % Combining the PDE bases should not affect order of BCs
            eq_num = eq_num_1;
        end

        % Finally, add the terms to the already specified equation.
        PDE_out.(obj_1){eq_num}.term = [PDE_1_out.(obj_1){eq_num}.term, PDE_2_out.free{ii}.term];
        % If PDE_out.(obj_1) previously specified a zero equation, it no
        % longer does now
        if isfield(PDE_out.(obj_1){eq_num},'is_zero')                       % DJ, 01/03/2025;
            PDE_out.(obj_1){eq_num}.is_zero = false;
        end
    end

end

end