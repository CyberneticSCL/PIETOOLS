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
% DJ, 06/23/2024: Initial coding;
% DJ, 01/09/2025: Add support for temporal differentiation of vectors of
%                   state variables;
% DJ, 02/18/2026: Add support for nonlinear terms;

% % % Process the inputs

% % Check that the PDE structure is properly specified.
% Check that the PDE is specified as 'pde_struct'.
if ~isa(PDE_in,'pde_struct')
    error("The input PDE structure should correspond to a free term to be used in a PDE.")
end
% Make sure the PDE corresponds to a single state or set of terms.
[is_pde_var_in,obj] = is_pde_var(PDE_in);
if ~is_pde_term(PDE_in)
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
    elseif isscalar(dval)
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
        % Check that the input corresponds to a column of PDE state
        % variables
        for ii=1:numel(PDE_out.free)                                        % DJ, 01/09/2025
            if isfield(PDE_out.free{ii},'term') && ~isempty(PDE_out.free{ii}.term)
                if numel(PDE_out.free{ii}.term)>1
                    error('Temporal differentiation of multiple PDE terms is not supported');
                end
                trm = PDE_out.free{ii}.term{1};
                if numel(PDE_out.free{ii}.term{1})>1
                    error("Temporal differentiation of nonlinear terms is not supported.")
                end
                if ~isfield(trm,'x')
                    % Term must involve a state variable
                    error('Temporal differentiation of input or output signals is not supported.');
                % elseif isfield(trm,'C') && (~isdouble(trm.C) ||...
                %         (~isequal(double(trm.C),1) && ~isequal(double(trm.C),eye(PDE_in.x{1}.size))))
                %     % Variable cannot be multiplied with anything.
                %     error('Temporal differentiation of states multiplied with coefficients is not supported.');
                elseif isfield(trm,'D') && any(trm.D)
                    % Variable cannot be differentiated.
                    error('Simultaneous differentiation with respect to time and space is not supported.');
                elseif isfield(trm,'loc') && ~isempty(trm.loc) ...
                        && (any(size(trm.loc)~=size(PDE_in.x{trm.x}.vars'))...
                            || any(~isequal(polynomial(trm.loc),polynomial(PDE_in.x{trm.x}.vars'))))
                    % Variable cannot be evaluated at a particular position.
                    error('Evaluating temporal derivatives of states at a spatial position is not supported.');
                elseif isfield(trm,'delay') && (~isdouble(trm.delay) || ~double(trm.delay)==0)
                    % Variable cannot be delayed in time.
                    error('Temporal differentiation of delayed states is not supported.');
                elseif isfield(trm,'I') 
                    % Variable cannot be integrated.
                    for jj=1:numel(trm.I)
                        if ~isempty(trm.I{jj})
                            error('Temporal differentiation of integrals of states is not supported.');
                        end
                    end
                end
            end
        end
    end

    % Initialize a new term with temporal derivative of state variable
    % increased by dval
    for ii=1:numel(PDE_out.free)
        if isfield(PDE_out.free{ii},'term') && ~isempty(PDE_out.free{ii}.term)
            if isfield(PDE_out.free{ii}.term{1},'tdiff')
                PDE_out.free{ii}.term{1}.tdiff = PDE_out.free{ii}.term{1}.tdiff + dval;
            else
                PDE_out.free{ii}.term{1}.tdiff = dval;
            end
        end
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
        % If not, differentiating returns 0.
        PDE_out.free{1}.vars = polynomial(zeros(0,1));
        PDE_out.free{1}.size = PDE_in.([obj,'_tab'])(1,2);
        PDE_out.free{1}.term = cell(1,0);
        return;
        %error("Differentiation is only supported with respect to variables on which the state depends.")
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
    PDE_out.free{ii}.term = {};
    n_trms = numel(PDE_in.free{ii}.term);
    % % Loop over all the terms, taking the desired derivative.
    for jj=1:n_trms
        % Make sure no derivative is taken of an output, or the temporal
        % derivative of a state.
        term_jj = PDE_in.free{ii}.term{jj};
        if is_LHS_term(term_jj)
            error("Differentiation of an output object or temporal derivative of a state is not supported.")
        end
        % Make sure the term has a coefficient
        if ~isfield(term_jj,'C')
            term_jj(1).C = 1;
        end
        % In nonlinear case, we have to use the product rule.
        % Compute all orders of derivatives of the factors with respect to
        % each of the variables that add up to specified order
        n_fctrs = numel(term_jj);
        dval_fctr = zeros(1,0);
        nreps_fctr = 1;
        for var_num=1:nvars
            d = dval(var_num);
            dval_tmp = (0:d)';
            for j=2:n_fctrs
                dval_tmp = [repmat(dval_tmp,d+1,1),repelem((0:d)',size(dval_tmp,1),1)];
                dval_tmp = dval_tmp(sum(dval_tmp,2)<=d,:);
            end
            dval_tmp = dval_tmp(sum(dval_tmp,2)==d,:);            
            dval_fctr = [repmat(dval_fctr,size(dval_tmp,1),1),repelem(dval_tmp,size(dval_fctr,1),1)];
            % Determine how many times each product of terms appears
            % (based on multinomial theorem)
            nreps_tmp = factorial(d)./prod(factorial(dval_tmp),2);
            nreps_fctr = repmat(nreps_fctr,size(nreps_tmp,1),1).*repelem(nreps_tmp,size(nreps_fctr,1),1);
        end
        new_terms = cell(1,0);
        trm_num = 1;
        fctr_strt = ones(size(dval_fctr,1),1);
        d_idx = 1;
        vars_ii = [];
        while d_idx<=size(dval_fctr,1)
            % Loop over each of the obtained derivatives of the different
            % factors
            rm_term = false;
            if d_idx>numel(new_terms)
                new_terms = [new_terms,{term_jj}];
            end
            % Multiply the derivative with the number of times it appears
            new_terms{d_idx}(1).C = nreps_fctr(d_idx)*new_terms{d_idx}(1).C;
            % Loop over the different factors in the term, taking the
            % derivative of order dval_fctr(k,l) of the lth factor
            term_k = new_terms{d_idx};
            vars_k = [];
            for fctr_num = fctr_strt(d_idx):n_fctrs
                % Extract the factor and order of derivative to be taken
                fctr = term_k(fctr_num);
                dval_jj = dval_fctr(d_idx,fctr_num:n_fctrs:end);
                if all(dval_jj==0)
                    continue
                end
                % Make sure derivative is only taken of state variables, not inputs.
                if (isfield(fctr,'u') && ~isempty(fctr.u)) ||...
                    (isfield(fctr,'w') && ~isempty(fctr.w))
                    error('Differentiation of input signals is currently not supported.')
                end
                % Determine which state variable is involved in the factor.
                Ridx = fctr.x;
                % Determine which variables the state depends on.
                vars_obj = PDE_in.x{Ridx}.vars;
        
                % Check if the coefficients depend on the variables with which to
                % differentiate.
                Cjj = polynomial(fctr.C);
                if ~any(ismember(vars.varname,Cjj.varname))
                    % The coefficients need not be differentiated
                    % --> just take derivative of the state.
                    if any(~ismember(vars.varname,vars_obj.varname))
                        % The state does not depend on one of the variables
                        % --> differentiation will return zero.
                        rm_term = true;
                        break
                    else
                        % Otherwise, add order of derivative to current order of
                        % derivative.
                        dval_obj = compute_d_obj(vars_obj,vars,dval_jj);
                        new_terms{trm_num}(fctr_num) = fctr;
                        new_terms{trm_num}(fctr_num).D = fctr.D +dval_obj;
                        vars_k = unique([vars_k; vars_obj.varname]);
                    end
                else
                    % The coefficients need to be differentiated as well...
                    % --> use product rule.
        
                    % First, get rid of variables on which the state does not
                    % depend:
                    vars_jj = vars;
                    if any(~ismember(vars.varname,vars_obj.varname))
                        % If the state does not depend on certain variables,
                        % but C does, only the derivative of C wrt this 
                        % variable will be nonzero in the product rule.
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
                        if all(isequal(Cjj,0))
                            % Don't add zero term.
                            rm_term = true;
                            break
                        else
                            % Replace coefficients with their derivative.
                            fctr.C = Cjj;
                            % Keep track of which variables the term depends on.
                            tmp_vars_ii_jj = addvars(vars_obj,Cjj);     
                            vars_k = unique([vars_k;tmp_vars_ii_jj.varname]);
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
                    new_terms{trm_num}(fctr_num) = fctr;
                    new_terms{trm_num}(fctr_num).D = fctr.D +dval_obj;
        
                    % Then, add all other terms D^(i)C *D^(dval-i)x            
                    n_prods = prod(dval_obj+1);
                    trm_num_ll = trm_num+1;
                    for ll=2:n_prods
                        % Determine which order derivative is taken of C
                        dval_C = cell(1,length(dval_obj)+1);
                        [dval_C{:}] = ind2sub([dval_obj+1,1],ll);
                        dval_C = cell2mat(dval_C(1:end-1)) - 1;
        
                        % Differentiate coefficients with the variable.
                        Cll = polydiff(Cjj,vars_obj,dval_C);
                        if all(isequal(Cll,0))
                            % Don't add a zero term to the PDE...
                            continue
                        end
                        % Multiply to match how many times this combination of
                        % derivative of C and of x appears.
                        Cfctr = prod(factorial(dval_obj)./(factorial(dval_obj-dval_C).*factorial(dval_C)));
                        Cll = Cfctr*Cll;
        
                        % Set a new term with the differentiated coefficients, and 
                        % lower-order derivative of state.
                        new_terms{trm_num_ll} = new_terms{d_idx};
                        new_terms{trm_num_ll}(fctr_num) = fctr;
                        new_terms{trm_num_ll}(fctr_num).D = fctr.D +(dval_obj-dval_C);
                        new_terms{trm_num_ll}(fctr_num).C = Cll;
        
                        % Keep track of which variables the term depends on
                        tmp_vars_ii_jj = addvars(vars_obj,Cll);     
                        vars_k = unique([vars_k;tmp_vars_ii_jj.varname]);
                    
                        % Proceed to the next term.
                        dval_fctr = [dval_fctr(1:d_idx,:); dval_fctr(d_idx,:); dval_fctr(d_idx+1:end,:)];
                        nreps_fctr = [nreps_fctr(1:d_idx); 1; nreps_fctr(d_idx+1:end)];
                        fctr_strt = [fctr_strt(1:d_idx); fctr_num+1; fctr_strt(d_idx+1:end)];
                        trm_num_ll = trm_num_ll+1;
                    end
                    
                end
            end
            if rm_term
                % The current derivative produces 0, move on to the next
                dval_fctr(d_idx,:) = [];
                nreps_fctr(d_idx) = [];
                new_terms(d_idx) = [];
                continue
            else
                % Move on to the next derivative
                d_idx = d_idx+1;
                trm_num = trm_num+1;
                vars_ii = unique([vars_ii; vars_k]);
            end
        end
        PDE_out.free{ii}.term = [PDE_out.free{ii}.term, new_terms];
    end
    % Keep track of which variables the equation depends on.
    PDE_out.free{ii}.vars = polynomial(vars_ii);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function dval_obj = compute_d_obj(vars_obj,vars,dval)
% Compute the order of the derivatives to be taken of variables "vars_obj",
% given the desired order of the derivatives "dval" of the variables
% "vars".

dval_obj = zeros(1,size(vars_obj,1));
for jj=1:size(vars_obj,1)
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