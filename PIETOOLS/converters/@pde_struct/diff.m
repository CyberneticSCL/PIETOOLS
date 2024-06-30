function PDE_out = diff(PDE_in,vars,dval)
% PDE_OUT = DIFF(PDE_IN,VARS,D) declares a new PDE term corresponding
% the derivative of the input PDE term "PDE_in" with respect to variables
% "vars" and up to order "dval".
%
% INPUT
% - PDE_in:     'pde_struct' object corresponding to either a single PDE
%               variable (state, input, or output) or corresponding
%               to one or several free terms in a PDE. Differentiation is
%               only supported for state variables!
% - vars:       nx1 array of type 'polynomial' with each element speicfying
%               a single variable (pvar object) with which to
%               differentiate.
%               Can also be set to 't' (char object) or the pvar 't' to
%               take a temporal derivative of the state variable.
% - dval:       nx1 array of type 'double' with each element a nonnegative
%               integer specifying the desired order of the derivative of
%               the terms in the PDE structure with respect to the
%               associated variable in "vars".
%
% OUTPUT
% - PDE_out:    'pde_struct' object representing the same terms as the
%               input PDE structure, but differentiated up to the desired
%               order with respect to the desired variables.
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
elseif is_pde_var_in && ~strcmp(obj,'x')
    error("Differentiation of inputs or outputs is not supported.")
end

% % Check that the variables are properly specified.
if ischar(vars)
    vars = polynomial({vars});
elseif iscellstr(vars)
    vars = polynomial(vars);
elseif ~isa(vars,'polynomial') || ~ispvar(vars)
    error('The variables to differentiate with should be specified as nx1 array of pvar objects')
end
% Check the dimensions of the second input
if all(size(vars)>1)
    error('The variables with which to differentiate should be specified as nx1 array of pvar objects')
else
    % Convert to column array
    vars = vars(:);
    nvars = length(vars);
end
% Make sure the variables are unique
if length(vars.varname)~=nvars
    error('The variables with which to differentiate should be unique.')
end

% % Check that the orders of the derivatives are properly specified.
if nargin<=2
    % Assume first-order derivative is desired.
    dval = ones(nvars,1);
else
    if ~isnumeric(dval) || any(dval<0) || any(dval~=round(dval))
        error('The order of the derivative should be specified as nx1 array of nonnegative integer values.')
    end
    % Check the dimensions of the third input
    if all(size(dval)>1)
        error('The order of the derivative should be specified as nx1 array of nonnegative integer values.')
    elseif numel(dval)==1
        % Assume same order derivative with respect to all variables.
        dval = dval*ones(nvars,1);
    elseif length(dval)~=nvars
        error('The number of orders of derivatives should match the number of variables with which to differentiate.')
    else
        % Convert to column array.
        dval = dval(:);
    end
    % Get rid of 0-order derivatives.
    vars = vars(dval~=0);       dval = dval(dval~=0);
    nvars = length(vars);
end

% Initialize the output structure;
PDE_out = PDE_in;


% % % Deal with the case of temporal differentiation.
if ismember('t',vars.varname)
    if nvars>1
        error("Differentiation with respect to temporal and spatial variable simultaneously is not supported.")
    elseif ~is_pde_var_in
        error('Differentiation of PDE terms with respect to temporal variable is not supported');
    end
    % Initialize a new term with temporal derivative of state variable
    % equal to dval
    for ii=1:numel(PDE_out.x)
        PDE_out.free{1}.size = PDE_in.x{ii}.size;
        PDE_out.free{1}.vars = PDE_in.x{ii}.vars;

        PDE_out.free{1}.term{1}.(obj) = ii;      % Index of state variable in the structure, not its ID!
        PDE_out.free{1}.term{1}.tdiff = dval;
        PDE_out.free{1}.term{1}.C = 1;
    end
    return
end

% % % Deal with the case that the input PDE corresponds to just a state 
% % % variable
% % % --> Initialize a PDE term with the desired derivative.
if is_pde_var_in
    % Extract spatial variables on which the PDE object depends.
    vars_obj = PDE_in.(obj){1}.vars;    nvars_obj = size(vars_obj,1);

    % Check that the object indeed depends on the variables with which to
    % differentiate.
    if any(~ismember(vars.varname,vars_obj.varname))
        error("Differentiation is only supported with respect to variables on which the state depends.")
    end

    % Initialize a PDE term with no derivative
    PDE_out.free{1}.term{1}.(obj) = 1; % or should it be PDE_in.([obj,'_tab'])(1,1);
    PDE_out.free{1}.term{1}.C = eye(PDE_in.([obj,'_tab'])(1,2));
    PDE_out.free{1}.term{1}.loc = vars_obj';
    PDE_out.free{1}.term{1}.I = cell(nvars_obj,1);

    % Set the values of the derivative.
    PDE_out.free{1}.term{1}.D = compute_d_obj(vars_obj,vars,dval);

    % Keep track of which variables the newlt started equation depends on.
    PDE_out.free{1}.vars = vars_obj;

    % Keep track of the size of the equation we are building.
    PDE_out.free{1}.size = PDE_in.([obj,'_tab'])(1,2);

    return
end

% % % Finally, proceed with case where derivative is taken of one or
% % % multiple terms.
% % Loop over all equations, and take the derivative.
for ii=1:numel(PDE_out.free)
    PDE_out.free{ii}.trm = {};
    n_trms = numel(PDE_in.free{ii}.term);
    trm_num = 1;
    % % Loop over all the terms, taking the desired derivative.
    for jj=1:n_trms
        % Make sure no derivative is taken of an output, or the temporal
        % derivative of a state.
        if is_LHS_term(PDE_out.free{ii}.term{jj})
            error("Differentiation of an output object or temporal derivative of a state is not supported.")
        end
        % Make sure derivative is only take of state variables, not inputs.
        if isfield(PDE_out.free{ii}.term{jj},'x')
            Robj = 'x';
        else
            error('Differentiation of input signals is currently not supported.')
        end
        % Determine which state variable is involved in the term.
        %Ridx = PDE_out.([obj,'_tab'])(:,1)==PDE_out.free{ii}.term{jj}.obj;
        Ridx = PDE_out.free{ii}.term{jj}.(obj);
        % Determine which variables the state depends on.
        vars_obj = PDE_out.(obj){Ridx}.vars;

        % Check if the coefficients depend on the variables with which to
        % differentiate.
        Cjj = polynomial(PDE_in.free{ii}.term{jj}.C);
        if ~any(ismember(vars.varname,Cjj.varname))
            % The coefficients need not be differentiated
            % --> just take derivative of the state.
            if any(~ismember(vars.varname,vars_obj.varname))
                % The state does not depend on one of the variables
                % --> differentiation will return zero.
                continue
            else
                % Otherwise, add order of derivative to current order of
                % derivative.
                dval_obj = compute_d_obj(vars_obj,vars,dval);
                PDE_out.free{ii}.term{trm_num} = PDE_in.free{ii}.term(jj);
                PDE_out.free{ii}.term{trm_num}.D = PDE_in.free{ii}.term(jj).D +dval_obj;
                trm_num = trm_num+1;
            end
        else
            % The coefficients need to be differentiated as well...
            % --> use product rule.

            % First, get rid of variables on which the state does not
            % depend:
            vars_jj = vars;     dval_jj = dval;
            tmp_vars_ii = [];
            if any(~ismember(vars.varname,vars_obj.varname))
                % If the state does not depend on certain variables, but C
                % does, only the derivative of C wrt this variable will be
                % nonzero in the product rule.
                for kk=1:nvars
                    if ~any(isequal(vars(kk),vars_obj))
                        % The object does not depend on variable kk
                        for ll=1:dval_jj(kk)
                            % Have to perform differentiation of
                            % 'polynomial' object one order at a time...
                            Cjj = diff(Cjj,vars_jj(kk));
                        end
                        dval_jj(kk) = 0;
                    end
                end
                if all(isequal(PDE_out.free{ii}.term{trm_num}.C,0))
                    % Don't add zero term.
                    continue
                else
                    % Add remaining nonzero term from product rule.
                    PDE_out.free{ii}.term{trm_num} = PDE_in.free{ii}.term{jj};
                    PDE_out.free{ii}.term{trm_num}.C = Cjj;
                    % Keep track of which variables the term depends on.
                    tmp_vars_ii_jj = addvars(vars_obj,Cjj);     
                    tmp_vars_ii = unique([tmp_vars_ii;tmp_vars_ii_jj.varname]);
                end
                % Get rid of zero-order derivatives.
                vars_jj = vars_jj(dval_jj~=0);
                dval_jj = dval_jj(dval_jj~=0);
            end
            
            % % For the rest, we need a repeated product rule:
            % (d/ds1)^k1 ... (d/dsn)^kn (C*x)
            % = sum_{i1=0}^k1 nchoosek(k1,i1)
            %    * sum_{i2=0}^{k2} nchoosek(k2,i2)
            %      ...
            %       * sum_{in=0}^{kn} nchoose(kn,in)
            %         * (d/ds1)^i1 ... (d/dsn)^in C
            %         * (d/ds1)^(k1-i1) ... (d/dsn)^(kn-in) x
            % % Add each of these terms to the PDE

            % First account for term C*D^(dval)x
            dval_obj = compute_d_obj(vars_obj,vars_jj,dval_jj);
            PDE_out.free{ii}.term{trm_num}.D = PDE_in.free{ii}.term{jj}.D +dval_obj;
            trm_num = trm_num+1;

            % Then, add all other terms D^(i)C *D^(dval-i)x            
            n_prods = prod(dval_obj+1);
            for ll=2:n_prods
                % Determine which order derivative is taken of C
                dval_C = cell(1,length(dval_obj));
                [dval_C{:}] = ind2sub(dval_obj+1,ll);
                dval_C = cell2mat(dval_C) - 1;

                % Differentiate coefficients with the variable.
                Cll = polydiff(Cjj,vars_obj,dval_C);
                if all(isequal(Cll,0))
                    % Don't add a zero term to the PDE...
                    continue
                end
                % Multiply to match how many times this combination of
                % derivative of C and of x appears.
                fctr = prod(factorial(dval_obj)./(factorial(dval_obj-dval_C).*factorial(dval_C)));
                Cll = fctr*Cll;

                % Set a new term with the differentiated coefficients, and 
                % lower-order derivative of state.
                PDE_out.free{ii}.term{trm_num} = PDE_in.free{ii}.term{jj};
                PDE_out.free{ii}.term{trm_num}.D = PDE_in.free{ii}.term{jj}.D +(dval_obj-dval_C);
                PDE_out.free{ii}.term{trm_num}.C = Cll;

                % Keep track of which variables the term depends on
                tmp_vars_ii_jj = addvars(vars_obj,Cll);     
                tmp_vars_ii = unique([tmp_vars_ii;tmp_vars_ii_jj.varname]);
            
                % Proceed to the next term.
                trm_num = trm_num+1;
            end
            % Keep track of which variables the equation depends on.
            PDE_out.free{ii}.vars = polynomial(tmp_vars_ii);
        end
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function dval_obj = compute_d_obj(vars_obj,vars,dval)
% Compute the order of the derivatives to be taken of variables "vars_obj",
% given the desired order of the derivatives "dval" of the variables
% "vars".

dval_obj = zeros(1,length(vars_obj));
for jj=1:length(vars_obj)
    % Find variable "vars_obj(jj)" in "vars".
    idx = isequal(vars_obj(jj),vars);
    % Set the associated order of differentiability.
    if any(idx)
        dval_obj(jj) = dval(idx);
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function newvars = addvars(oldvars,C)
% Add variables that appear in polynomial object C to the list of variables
% oldvars.

if isdouble(C)
    % C introduces no new variables.
    newvars = oldvars;
else
    % Combine spatial variabels of the PDE variable and the factor C
    % Has to be done somewhat akwardly, as the polynomial structure
    % rearranges the variables alphabetically
    nvars = size(oldvars,1);
    vnames = cell(nvars+length(C.varname),1);
    for jj=1:nvars
        vnames{jj} = oldvars(jj,1).varname{1};
    end
    vnames(nvars+1:end) = C.varname;
    vnames = unique(vnames,'stable');
    newvars = polynomial(vnames);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function dC = polydiff(C,vars,dval)
% Compute derivative of 'polynomial' class object "C" with each of the
% variables "vars(jj)" up to order "dval(jj)";

% Very basic (expensive) implementation, just use for loops
dC = C;
for ii=1:length(vars)
    for jj=1:dval(ii)
        dC = diff(dC,vars(ii));
    end
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %