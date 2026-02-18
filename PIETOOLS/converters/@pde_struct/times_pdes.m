function PDE_out = times_pdes(PDE1,PDE2)
% PDE_OUT = TIMES_PDES(PDE1,PDE2) computes the elementwise product of two 
% terms that may appear in a PDE.
%
% INPUTS
% - PDE1:   'pde_struct' object representing a single term or set of terms
%           to appear in a PDE;
% - PDE2:   'pde_struct' object representing a single term or set of terms
%           to appear in a PDE. The number of rows of terms must match that
%           in PDE1, and the vector-dimension of each row of terms must
%           match that in PDE1 as well;
%
% OUTPUTS
% - PDE_out:    'pde_struct' object representing the elementwise product of
%               the sets of input terms.

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
% DJ, 02/13/2026: Initial coding;

% Make sure the inputs are 'pde_struct' objects
if ~isa(PDE1,'pde_struct') || ~isa(PDE2,'pde_struct')
    error("The factors to multiply should be specified as 'pde_struct' objects.")
end
% Make sure the input corresponds to a single state/input/output
% or set of terms
[is_pde_var1,obj1] = is_pde_var(PDE1);
[is_pde_var2,obj2] = is_pde_var(PDE2);
if (~is_pde_var1 && ~is_pde_term(PDE1)) || (~is_pde_var2 && ~is_pde_term(PDE2))
    error("The factors to multiply should correspond to individual terms in the PDE.")
end
% Convert the PDE variables to terms
if is_pde_var1
    PDE1 = var2term(PDE1,obj1);
end
if is_pde_var2
    PDE2 = var2term(PDE2,obj2);
end
% Express the terms in terms of the same set of states, inputs, and
% outputs
[PDE1,PDE2] = pde_common_basis(PDE1,PDE2);

% Extract the number of rows in the different terms
eq_size_list1 = size(PDE1,'free','vec_size');
eq_size_list2 = size(PDE2,'free','vec_size');
% Augment scalar term to vector for elementwise multiplication
if isscalar(eq_size_list1) && ~isscalar(eq_size_list2)
    PDE1.free = repmat(PDE1.free,[numel(eq_size_list2),1]);
    eq_size_list1 = repmat(eq_size_list1,[numel(eq_size_list2),1]);
elseif isscalar(eq_size_list2) && ~isscalar(eq_size_list1)
    PDE2.free = repmat(PDE2.free,[numel(eq_size_list1),1]);
    eq_size_list2 = repmat(eq_size_list2,[numel(eq_size_list1),1]);
end
% Make sure the vector dimensions match for multiplication
if numel(eq_size_list1)~=numel(eq_size_list2)
    error("Number of rows of terms must match for elementwise multiplication.")
end
is_neq = eq_size_list1~=eq_size_list2;
if any(eq_size_list1(is_neq)>1 & eq_size_list2(is_neq)>1)
    error("Row dimensions of terms must match for elementwise multiplication.")
end
eq_size_out = max(eq_size_list1,eq_size_list2);

% Construct a structure representing the product
PDE_out = PDE1;
for i=1:numel(PDE_out.free)
    % Set the ith row of terms in the PDE
    eq_i = PDE_out.free{i};
    eq1_i = PDE1.free{i};
    eq2_i = PDE2.free{i};
    % The product will depend on the variables of both factors
    vars1 = eq1_i.vars;
    vars2 = eq2_i.vars;
    new_varnames = unique([vars1.varname;vars2.varname]);
    eq_i.vars = polynomial(new_varnames);
    % The size of the elementwise product will be the size of the factors
    eq_i.size = eq_size_out(i);
    % The terms in the ith row will be the products of the terms in 1 and 2
    eq_i.term = cell(1,numel(eq1_i.term)*numel(eq2_i.term));
    for j=1:numel(eq_i.term)
        [trm_num1,trm_num2] = ind2sub([numel(eq1_i.term),numel(eq2_i.term)],j);
        % Do not allow for multiplication of terms on the left-hand side of the PDE
        if is_LHS_term(eq1_i.term{trm_num1}) || is_LHS_term(eq2_i.term{trm_num2})
            error("Multiplication of terms appearing on the left-hand side of the PDE is not supported.")
        end
        % To represent the product of two terms, we concatentate the
        % struct objects        
        eq_i.term{j} = times_term(eq1_i.term{trm_num1},eq2_i.term{trm_num2},PDE_out);
    end
    PDE_out.free{i} = eq_i;
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function term_out = times_term(term1,term2,PDE)
% TERM_OUT = TIMES_TERM(TERM1,TERM2,PDE) computes the product of two terms
% in a PDE.
%
% INPUTS
% - term1:  1xm struct representing a single term in a PDE as per the
%           'pde_struct' format;
% - term2:  1xm struct representing a single term in a PDE as per the
%           'pde_struct' format;
% - PDE:    'pde_struct' object including the information on the different
%           state, input, and output variables that may appear in the two
%           terms to multiply;
%
% OUTPUTS
% - term_out:   1xm*n struct representing the product of the specified
%               terms. This struct will be given by the concatenation of
%               the structs term1 and term2;
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
% DJ, 02/13/2026: Initial coding;

% Check whether term1 involves a state or input variable
if isfield(term1,'x')
    obj1 = 'x';
else
    error("Multiplication with input or output signals is not supported.")
end
obj1_num = term1.(obj1);
obj1_size = PDE.(obj1){obj1_num}.size;
obj1_vars = PDE.(obj1){obj1_num}.vars;

% Check whether term2 involves a state or input variable
if isfield(term2,'x')
    obj2 = 'x';
else
    error("Multiplication with input or output signals is not supported.")
end
obj2_num = term2.(obj2);
obj2_size = PDE.(obj1){obj2_num}.size;
obj2_vars = PDE.(obj2){obj2_num}.vars;

% Add potential missing object field to each term (necessary for
% concatenation)
if ~strcmp(obj1,obj2)
    term1.(obj2) = [];
    term2.(obj1) = [];
end

% Preclude delay in the terms
if isfield(term1,'delay')
    if ~isequal(term1.delay,0)
        error("Multiplication of delayed signals is not supported.")
    else
        term1 = rmfield(term1,'delay');
    end
end
if isfield(term2,'delay')
    if ~isequal(term2.delay,0)
        error("Multiplication of delayed signals is not supported.")
    else
        term2 = rmfield(term2,'delay');
    end
end

% Check if the number of rows in the terms matches
if isfield(term1,'C')
    obj1_size = 1;
    if ~isscalar(term1(1).C)
        obj1_size = max(obj1_size,size(term1(1).C,1));
    end
else
    term1.C = eye(obj1_size);
end
if isfield(term2,'C')
    if ~isscalar(term2(1).C)
        obj2_size = size(term2(1).C,1);
    end
else
    term2.C = eye(obj2_size);
end
% Augment scalar value to vector
if obj2_size>1 && obj1_size==1
    term1.C = ones(obj2_size,1)*term1.C;
elseif obj1_size>1 && obj2_size==1
    term2.C = ones(obj1_size,1)*term2.C;
elseif obj1_size~=obj2_size
    error("Row dimensions of terms must match for elementwise multiplication.")
end

% Make sure the derivative to be taken is specified
if ~isfield(term1,'D')
    term1.D = zeros(1,size(obj1_vars,1));
end
if ~isfield(term2,'D')
    term2.D = zeros(1,size(obj2_vars,1));
end
% Make sure the location at which to evaluate is specified
if ~isfield(term1,'loc')
    term1.loc = obj1_vars(:,1)';
end
if ~isfield(term2,'loc')
    term2.loc = obj2_vars(:,1)';
end
% Make sure the integral to take is specified
if ~isfield(term1,'I')
    term1.I = cell(size(obj1_vars,1),1);
end
if ~isfield(term2,'I')
    term2.I = cell(size(obj2_vars,1),1);
end

% Finally, just concatenate the two structures
term_out = [term1,term2];

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %