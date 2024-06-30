function PDE_out = int(PDE_in,vars,L,U)
% PDE_OUT = INT(PDE_IN,VARS,L,U) declares a new PDE term corresponding
% the input PDE term "PDE_in" integrated from "L" up to evaluated "U" with
% respect to variables "vars".
%
% INPUT
% - PDE_in:     'pde_struct' object corresponding to either a single PDE
%               variable (state, input, or output) or corresponding
%               to one or several free terms in a PDE.
% - vars:       nx1 array of type 'polynomial' with each element speicfying
%               a single variable (pvar object) with which to integrate.
% - L:          nx1 array of type 'double' or 'polynomial', specifying the
%               lower limits of the integral for each variable in "vars".
%               Can also be nx2 array, in which case L(:,2) will be taken
%               as input U.
% - U:          nx1 array of type 'double' or 'polynomial', specifying the
%               upper limits of the integral for each variable in "vars".
%
% OUTPUT
% - PDE_out:    'pde_struct' object representing the same terms as the
%               input PDE structure, but with integration performed with
%               respect to variables "vars", from lower limits "L up to
%               upper limits "U".
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
elseif is_pde_var_in && ~strcmp(obj,'x') && ~strcmp(obj,'u') && ~strcmp(obj,'w')
    error("Integration of outputs is currently not supported.")
end

% % Check that the variables are properly specified.
if ischar(vars)
    vars = polynomial({vars});
elseif iscellstr(vars)
    vars = polynomial(vars);
elseif ~isa(vars,'polynomial') || ~ispvar(vars)
    error('The variables with which to integrate should be specified as nx1 array of pvar objects')
end
% Check the dimensions of the second input
if all(size(vars)>1)
    error('The variables with which to integrate should be specified as nx1 array of pvar objects')
else
    % Convert to column array
    vars = vars(:);
    nvars = length(vars);
end
% Make sure the variables are unique
if length(vars.varname)~=nvars
    error('The variables with which to integrate should be unique.')
end

% % Check that the location at which to evaluate is properly specified.
if nargin<=3 && size(L,2)==2
    % Allow domain of integration to be specified as [L,U].
    U = L(:,2);
    L = L(:,1);
elseif nargin<=3
    error("No upper limit has been specified for the integral.")
end
% Check the dimensions of the third input
if (~isa(L,'polynomial') && ~isa(L,'double')) || (~isa(U,'polynomial') && ~isa(U,'double'))
    error("The lower and upper limits of the integral should be specified as nx1 arrays of type 'polynomial' or 'double'.")
elseif all(size(L)>1) || all(size(U)>1)
    error("The lower and upper limits of the integral should be specified as nx1 arrays of type 'polynomial' or 'double'.")
elseif length(L)~=nvars || length(U)~=nvars
    error("The number of upper and lower limits of the integral should match the number of variables to integrate with.")
else
    % Convert to column array.
    L = L(:);       U = U(:);
end
% Make sure that each limit is either a position ('double') or variable
% ('pvar')
L = polynomial(L);      U = polynomial(U);
is_full_int = true(nvars,1);        % Keep track of for which variables a full integral is taken.
for kk=1:nvars
    if (~isdouble(L(kk)) && ~ispvar(L(kk))) || (~isdouble(U(kk)) && ~ispvar(U(kk)))
        error("Upper and lower limits of integrals should correspond to real values or single variables ('pvar' objects).")
    elseif ispvar(L(kk)) && ispvar(U(kk))
        error("Integration with both a variable upper and lower limit is not supported.")
    elseif isequal(L(kk),U(kk))
        error("Upper and lower limits of the integral should be distinct.")
    end
    if ispvar(L(kk)) || ispvar(U(kk))
        % If one of the limits is a variable, this is a partial integral.
        is_full_int = false;
    end
end
% Also check no integration is performed of temporal variable.
if ismember('t',vars.varname)
    error("Temporal delay is currently not supported.")
end

% % % Perform the actual integration.
% Convert single PDE variable (state, input) to a term.
if is_pde_var_in
    PDE_in = var2term(PDE_in);
end

% % Loop over all equations in the PDE, and take the integral of all terms.
PDE_out = PDE_in;
for ii=1:numel(PDE_out.free)
    if ~isfield(PDE_out.free{ii},'term') || isempty(PDE_out.free{ii}.term)
        % If there are no terms, there is nothing to be done.
        continue
    end
    for jj=1:numel(PDE_out.free{ii}.term)
        term_jj = PDE_out.free{ii}.term{jj};
        % First, check the term does not correspond to a left-hand side
        % object (a temporal derivative of state, or an output)
        if is_LHS_term(term_jj)
            error("Integration of outputs or temporal derivatives of state variables is not supported.")
        end
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
        % Initialize default integral if no location has been specified.
        obj_vars = PDE_in.(obj_jj){term_jj.(obj_jj)}.vars(:,1);
        if ~isfield(term_jj,'I')
            % Empty cell corresponds to no integration.
            term_jj.I = cell(length(obj_vars),1);
        end
        % Also check whether the object is already evaluated somewhere.
        if isfield(term_jj,'loc')
            term_vars = polynomial(term_jj.loc)';
        else
            % If no location is specified, evaluate at s=s.
            term_vars = obj_vars;
            term_jj.loc = obj_vars';
        end

        % Next, integrate for each variable separately.
        for ll=1:nvars
            % Check if the term actually depends on the variable
            is_term_var = isequal(vars(ll),term_vars);
            if any(is_term_var)
                % Check that the desired limits of the integral make sense.
                if (ispvar(L(ll)) && ~strcmp(L(ll).varname{1},obj_vars(is_term_var).varname{1})) ||...
                        (ispvar(U(ll)) && ~strcmp(U(ll).varname{1},obj_vars(is_term_var).varname{1}))
                    error("Variable limits of integral should match variable on which the component depends.")
                end
                % Set the desired integral limits.
                term_jj.I{is_term_var} = [L(ll),U(ll)];
                % If partial integration is performed, replace variable
                % with dummy variable.
                if ispvar(L(ll)) || ispvar(U(ll))
                    if strcmp(vars(ll).varname{1},obj_vars(is_term_var).varname{1})
                        % Generate a dummy variable.
                        varname_ll = vars(ll).varname{1};
                        dumvar = polynomial({[varname_ll,'_d']});
                        % Evaluate the term at the dummy variable.
                        term_jj.loc(is_term_var) = dumvar;
                        if isfield(term_jj,'C')
                            % Make sure to also replace variable in
                            % coefficients with dummy variable.
                            term_jj.C = subs(polynomial(term_jj.C),vars(ll),dumvar);
                        end
                    end
                end
            else
                % The term does not depend on the specified variable
                % --> only coefficients need to be integrated.
                if isfield(term_jj,'C')
                    term_jj.C = int(term_jj.C,vars(ll),L(ll),U(ll));
                else
                    term_jj.C = polynomial(U(ll)-L(ll));
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
            is_int_var = isequal(vars_ii(ll),vars);
            if any(is_int_var) && is_full_int(is_int_var)
                % If a full integral is taken, dependence on variable is
                % removed.
                retain_vars(ll) = false;
            end
        end
        PDE_out.free{ii}.vars = vars_ii(retain_vars);
    end
end

end