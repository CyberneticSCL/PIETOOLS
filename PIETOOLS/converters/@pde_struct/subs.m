function PDE_out = subs(PDE_in,vars,locs)
% PDE_OUT = SUBS(PDE_IN,VARS,LOC) declares a new PDE term corresponding
% the input PDE term "PDE_in" evaluated at vars=loc.
%
% INPUT
% - PDE_in:     'pde_struct' object corresponding to either a single PDE
%               variable (state, input, or output) or corresponding
%               to one or several free terms in a PDE.
% - vars:       nx1 array of type 'polynomial' with each element speicfying
%               a single variable (pvar object) which to evaluate at a
%               certain position.
% - locs:       nx1 array of type 'double' or 'polynomial' with each element 
%               specifying the desired real value or variable to substitute
%               the associated variable in "vars" with in the terms in the
%               PDE.
%
% OUTPUT
% - PDE_out:    'pde_struct' object representing the same terms as the
%               input PDE structure, but with the variables in "vars"
%               substituted by the values in "locs".
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
% Initial coding DJ - 06/23/2024


% % % Process the inputs

% % Check that the PDE structure is properly specified.
% Check that the PDE is specified as 'pde_struct'.
if ~isa(PDE_in,'pde_struct')
    error("The input PDE structure should correspond to a free term to be used in a PDE.")
end
% Make sure the PDE corresponds to a single state or set of terms.
[is_pde_var_in,obj] = is_pde_var(PDE_in);
if ~is_pde_var(PDE_in) && ~is_pde_term(PDE_in)
    error("The input PDE structure should correspond to a free term to be used in a PDE.")
elseif is_pde_var_in && ~(strcmp(obj,'x') || strcmp(obj,'u') || strcmp(obj,'w'))
    error("Substitution of outputs is currently not supported.")
end

% % Check that the variables are properly specified.
if ischar(vars)
    vars = polynomial({vars});
elseif iscellstr(vars)
    vars = polynomial(vars);
elseif ~isa(vars,'polynomial') || ~ispvar(vars)
    error('The variables to substitute should be specified as nx1 array of pvar objects')
end
% Check the dimensions of the second input
if all(size(vars)>1)
    error('The variables to substitute should be specified as nx1 array of pvar objects')
else
    % Convert to column array
    vars = vars(:);
    nvars = length(vars);
end
% Make sure the variables are unique
if length(vars.varname)~=nvars
    error('The variables to substitute should be unique.')
end

% % Check that the location at which to evaluate is properly specified.
% Check the dimensions of the third input
if ~isa(locs,'polynomial') && ~isa(locs,'double')
    error("The position at which to evaluate the terms should be specified as nx1 array of type 'polynomial' or 'double'.")
elseif all(size(locs)>1)
    error("The position at which to evaluate the terms should be specified as nx1 array of type 'polynomial' or 'double'.")
% elseif isa(locs,'double') && numel(locs)==1
%     % Assume same position for all variables.
%     locs = locs*ones(nvars,1);
elseif length(locs)~=nvars
    error("The number of positions at which to evaluate should match the number of variables to substitute.")
else
    % Convert to column array.
    locs = locs(:);
end
% Make sure that each element is either a position ('double') or variable
% ('pvar')
locs = polynomial(locs);
is_fixed_loc = true(nvars,1);       % Keep track of which variables are replaced with an actual position.
for kk=1:nvars
    if ~isdouble(locs(kk)) && ~ispvar(locs(kk))
        error("New objects with which to substitute the variables should correspond to real values or single variables ('pvar' objects).")
    end
    if ispvar(locs(kk))
        % Variable is replaced with a new variable.
        is_fixed_loc(kk) = false;
    end
end
% Also check no substitution is performed of temporal variable.
if ismember('t',vars.varname)
    error("Temporal delay is currently not supported.")
end


% % % Perform the actual substitution.
% Convert single PDE variable (state, input, output) to a single term.
if is_pde_var_in
    PDE_in = var2term(PDE_in);
end

PDE_out = PDE_in;
for ii=1:numel(PDE_out.free)
    if ~isfield(PDE_out.free{ii},'term') || isempty(PDE_out.free{ii}.term)
        continue
    end
    for jj=1:numel(PDE_out.free{ii}.term)
        term_jj = PDE_out.free{ii}.term{jj};
        % First, check the term does not correspond to a left-hand side
        % object (a temporal derivative of state, or an output)
        if is_LHS_term(term_jj)
            error("Substitution of outputs or temporal derivatives of state variables is not supported.")
        end
        % % Also check the term does not take an input, for which we also
        % % cannot perform substitution.
        % if isfield(term_jj,'u') || isfield(term_jj,'w')
        %     error("Substitution of inputs is not supported.")
        % end
        % Check what kind of object the term does correspond to.
        if isfield(term_jj,'x')
            obj_jj = 'x';
        elseif isfield(term_jj,'u')
            obj_jj = 'u';
        elseif isfield(term_jj,'w')
            obj_jj = 'w';
        else
            error("Integration of output signals is not supported.")
        end
        % Initialize default location if no location has been specified.
        obj_vars = PDE_in.(obj_jj){term_jj.(obj_jj)}.vars(:,1);
        if ~isfield(term_jj,'loc')
            term_jj.loc = obj_vars';
        end

        % Next, perform substitution for each variable separately.
        for ll=1:nvars
            % Substitute the variable in the coefficients.
            term_jj.C = subs(polynomial(term_jj.C),vars(ll),locs(ll));
            % Establish which of the variables in the state should be
            % substituted.
            isvar_ll = isequal(vars(ll),obj_vars);
            if ~any(isvar_ll)
                % The state does not depend on the variable
                % --> move on.
                continue
            else
                % The state does depend on the variable.
                % If the term is indeed evaluated at this variable,
                % replace.
                var_idx = find(isvar_ll);
                if isequal(vars(ll),term_jj.loc(var_idx))
                    term_jj.loc(var_idx) = locs(ll);
                end
                % If partial integral is performed, replace that variable
                % too.
                if isfield(term_jj,'I') && numel(term_jj.I)>=var_idx 
                    if ~isempty(term_jj.I{var_idx}) && ispolynomial(term_jj.I{var_idx})
                        term_jj.I{var_idx} = subs(term_jj.I{var_idx},vars(ll),locs(ll));
                    end
                end
            end
        end
        % Update the term in the PDE.
        PDE_out.free{ii}.term{jj} = term_jj;
    end
    % Update the list of variables on which the equation depends.
    if isfield(PDE_out.free{ii},'vars')
        vars_ii = PDE_out.free{ii}.vars;
        retain_vars = true(size(vars_ii,1),1);
        for ll=1:size(vars_ii,1)
            % Check if integration is performed with variable vars_ii(ll); 
            is_subs_var = isequal(vars_ii(ll),vars);
            if any(is_subs_var) && is_fixed_loc(is_subs_var)
                % If the variable is evaluated at a fixed position,
                % dependence on this variable is removed.
                retain_vars(ll) = false;
            end
        end
        PDE_out.free{ii}.vars = vars_ii(retain_vars);
    end
end

end