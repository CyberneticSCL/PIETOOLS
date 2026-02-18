function PDE_out = power(PDE_in,deg)
% PDE_OUT = POWER(PDE_IN,DEG) computes the elementwise power of a variable
% or term in a PDE.
%
% INPUTS
% - PDE_in: 'pde_struct' object representing a single term or set of m rows
%           of terms in a PDE;
% - deg:    m x 1 or scalar array specifying to which power to raise the
%           terms in each of the rows of the PDE;
%
% OUTPUTS
% - PDE_out:    'pde_struct' object representing the element-wise power of
%               the terms in the input PDE
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026  PIETOOLS Team
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
% DJ, 02/18/2026: Initial coding;


% Make sure the inputs are of appropriate type
if ~isa(PDE_in,'pde_struct')
    error("The first inputs must be a 'pde_struct' object.")
end
if ~isa(deg,'double') || any(deg<0) || any(deg~=round(deg))
    error("Degrees must be specified as real nonnegative integers.")
end
if any(deg==0)
    error("Degree-0 powers of PDE terms are currently not supported.")
end

% Make sure the input corresponds to a single state/input/output
% or set of terms
[is_pde_var_in,obj_in] = is_pde_var(PDE_in);
if (~is_pde_var_in && ~is_pde_term(PDE_in))
    error("Equations cannot be raised to a power.")
end
% Convert the PDE variable to a term
if is_pde_var_in
    PDE_in = var2term(PDE_in,obj_in);
end

% Extract the number of rows of terms
n_eqs = numel(PDE_in.free);
if isscalar(deg)
    deg = repmat(deg,[n_eqs,1]);
elseif numel(deg)~=n_eqs
    error("Number of degrees must match number of rows of terms.")
end

% Construct a structure representing the power
PDE_out = PDE_in;
for i=1:n_eqs
    % Set the ith row of terms in the PDE
    eq_i = PDE_out.free{i};
    n_terms = numel(eq_i.term);
    d = deg(i);
    if deg(i)==1 || n_terms==0
        % Raising to degree 1 returns the original terms
        continue
    end
    % Make sure each of the terms in the row is properly specified
    term_set = eq_i.term;
    for j=1:n_terms
        if is_LHS_term(term_set{j})
            error("Taking the power of terms appearing on the left-hand side of the PDE is not supported.")
        end
        % Fill in any missing fields in the struct defining this term   
        term_set{j} = complete_term(term_set{j},PDE_out);
    end
    % Compute all possible powers of the terms that add up to the desired
    % cumulative degree d
    degmat = (0:d)';
    for j=2:n_terms
        degmat = [repmat(degmat,d+1,1),repelem((0:d)',size(degmat,1),1)];
        degmat = degmat(sum(degmat,2)<=d,:);
    end
    degmat = degmat(sum(degmat,2)==d,:);
    % Determine with which factor to multiply each product of terms
    % (based on multinomial theorem)
    Cfctrs = factorial(d)./prod(factorial(degmat),2);
    % Take all products of terms that add up to degree d
    eq_i.term = cell(1,size(degmat,1));
    for k=1:size(degmat,1)
        % Concatenate d(k,j) copies of term j with the factor
        term_k = repmat(term_set{1},[1,degmat(k,1)]);
        for j=2:n_terms
            term_k = [term_k,repmat(term_set{j},[1,degmat(k,j)])];
        end
        % Multiply first factor with required number of repetitions
        term_k(1).C = Cfctrs(k)*term_k(1).C;
        eq_i.term{k} = term_k;
    end
    PDE_out.free{i} = eq_i;
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function term_out = complete_term(term_in,PDE)
% TERM_OUT = COMPLETE_TERM(TERM1,PDE) fills in the missing fields in a term
% in a PDE.
%
% INPUTS
% - term_in:    1xm struct representing a single term in a PDE as per the
%               'pde_struct' format;
% - PDE:    'pde_struct' object including the information on the different
%           state, input, and output variables that may appear in the term;
%
% OUTPUTS
% - term_out:   1xm struct representing the same term as the input, but now
%               with any missing fields filled in;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026  PIETOOLS Team
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
% DJ, 02/18/2026: Initial coding;

% If the term is a product of terms, we will assume all fields are filled
% in
if ~isscalar(term_in)
    term_out = term_in;
    return
end

% Make sure the term involves a state variable
if isfield(term_in,'x')
    obj = 'x';
else
    error("Taking power of input or output signals is not supported.")
end
obj_num = term_in.(obj);
obj_size = PDE.(obj){obj_num}.size;
obj_vars = PDE.(obj){obj_num}.vars;

% Preclude delay in the terms
term_out = term_in;
if isfield(term_out,'delay')
    if ~isequal(term_out.delay,0)
        error("Taking power of delayed signals is not supported.")
    else
        term_out = rmfield(term_out,'delay');
    end
end
% Make sure the coefficient matrix is specified
if ~isfield(term_out,'C')
    term_out.C = eye(obj_size);
end
% Make sure the derivative to be taken is specified
if ~isfield(term_out,'D')
    term_out.D = zeros(1,size(obj_vars,1));
end
% Make sure the location at which to evaluate is specified
if ~isfield(term_out,'loc')
    term_out.loc = obj_vars(:,1)';
end
% Make sure the integral to take is specified
if ~isfield(term_out,'I')
    term_out.I = cell(size(obj_vars,1),1);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %