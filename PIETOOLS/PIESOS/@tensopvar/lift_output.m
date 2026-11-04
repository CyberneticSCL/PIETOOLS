function Cop_out = lift_output(Cop,vars,dom,dims,deps)
% COP_OUT = LIFT_OUTPUT(COP,VARS,DOM,DIMS,DEPS) returns the 'tensopvar'
% object COP_OUT representing the tensor-PI operator COP with its output
% lifted to the function space defined by VARS.out, DOM.out and DEPS.out:
% the output is extended with variables it was not yet defined in, along
% which it is constant. For example, for Cop mapping to L2[a,b], and
% lifting to the space with variables (s,t) of degrees (1,1),
%   (Cop_out*x)(s,t) = (Cop*x)(s).
%
% INPUTS
% - Cop:    m x n 'tensopvar' object representing an operator
%               Cop: L2^n[dom_in] --> L2^m[dom_out];
% - vars:   1 x M cellstr specifying the names of the output variables of
%           the lifted operator. Must contain all output variables of Cop.
%           Can also be specified as struct with field 'out' containing
%           this cellstr (other fields, e.g. 'in', are ignored);
% - dom:    M x 2 'double' array specifying the intervals on which these
%           variables are defined. Must match the domains of the output
%           variables of Cop. Can also be specified as struct with field
%           'out' (other fields are ignored);
% - dims:   (optional) number of rows of the lifted operator, specified as
%           an integer, as 1 x 2 array [rows,columns], or as struct with
%           field 'out'. Defaults to the number of rows of Cop. If Cop is
%           scalar-valued (m=1), it can be lifted to a vector-valued
%           operator, repeating the scalar output in each row;
% - deps:   (optional) 1 x M array of nonnegative integers (or struct with
%           field 'out') specifying the degree of the lifted operator in
%           each of the output variables vars.out, i.e. the number of
%           copies of that variable in the output function space (see
%           Cop.deps). Must be at least the degrees of Cop, and defaults to
%           these degrees, or to 1 for a new variable;
%
% OUTPUTS
% - Cop_out: m x n 'tensopvar' object (or dims(1) x n, if specified)
%           representing the lifted operator, with
%               Cop_out.vars.out = sort(vars.out);
%               Cop_out.deps.out = deps.out (sorted accordingly);
%
% NOTES
% - The factors of the tensor-PI operator are lifted individually (by
%   'lift_output' for 'sopvar' objects). Each copy of a variable that is
%   added to the output is introduced in a factor that does not map to that
%   variable yet (Cop.depmat.out==0), taking the last such factors.
% - For a Kronecker product (Cop.otype(1)==2), the copies s_1,s_2,... of
%   an output variable are numbered over the factors in order (see
%   'apply_op'), so the factors that are lifted determine which copies
%   are constant. For example, for factors mapping to R and to L2[s],
%   lifting to deps=2 for s lifts the first factor, so that the output
%   is Y(s_2), with s_1 the constant copy.
% - For a pointwise product (Cop.otype(1)==1), the output can depend on a
%   new variable with degree 1 only (introduced in a single factor), and
%   the degrees in the existing variables cannot be increased.
% - Factors of type 'intop' (functionals) cannot be lifted. Factors that
%   are 'tensopvar' objects (nested operators) are lifted recursively.
% - A new variable with deps.out=0 is only added to the list of
%   variables (as in 'common_vars'), without changing the operator.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lift_output
%
% Copyright (C) 2026 PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 11/03/2026: Initial coding;

% Check that the tensor-PI operator is properly declared
if ~isa(Cop,'tensopvar')
    error("First argument must be of type 'tensopvar'.")
end
if ~Cop.chkval
    Cop = initialize(Cop);
end
if isempty(Cop.ops)
    error("Cannot lift an operator without terms.")
end
if nargin<=3
    dims = [];
end
[ntrms,d] = size(Cop.ops);
[vars_out,dom_out,m_new] = parse_inputs(vars,dom,dims,Cop.vars,Cop.dom,size(Cop));
M = numel(vars_out);

% Extract the dependency of the operator on the new output variables, with
% the degrees of the old variables, and 0 for the new variables
[~,idcs_old] = ismember(Cop.vars.out,vars_out);
deps_old = zeros(1,M);              deps_old(idcs_old) = Cop.deps.out;
depmat_old = zeros(d,M);            depmat_old(:,idcs_old) = Cop.depmat.out;

% Check that the degrees of the lifted operator are properly specified
if nargin<=4 || isempty(deps)
    deps = max(deps_old,1);
elseif isa(deps,'struct') && isfield(deps,'out')
    deps = deps.out;
end
if ~isnumeric(deps) || numel(deps)~=M || any(deps~=round(deps)) || any(deps<0)
    error("Degrees must be specified as 1 x M array of nonnegative integers, with M the number of output variables.")
end
deps = reshape(deps,1,[]);
if any(deps<deps_old)
    error("The degrees of the lifted operator cannot be smaller than those of the operator.")
end

% Establish in which factors to introduce the copies of the new variables,
% in the last factors that do not map to the variable yet. These must be
% factors that can be lifted
can_lift = true(d,1);
for j=1:d
    can_lift(j) = ~any(cellfun('isclass',Cop.ops(:,j),'intop'));
end
host = false(d,M);
for k=find(deps>deps_old)
    if Cop.otype(1)==1 && (deps_old(k)>0 || deps(k)>1)
        error("For a pointwise product, the output can only depend on a new variable with degree 1, and the degree in existing variables cannot be increased.")
    end
    idcs_j = find(can_lift & depmat_old(:,k)==0);
    nadd = deps(k)-deps_old(k);
    if numel(idcs_j)<nadd
        error("Cannot lift to degree %d in variable '%s': only %d factors can be lifted along this variable.",deps(k),vars_out{k},numel(idcs_j))
    end
    host(idcs_j(end-nadd+1:end),k) = true;
end

% Lift the factors
Cop_out = Cop;
for j=find(any(host,2))'
    is_k = host(j,:);
    for i=1:ntrms
        op_ij = Cop.ops{i,j};
        % Output variables of the factor, with the added variables
        vars_ij = struct('out',{[reshape(op_ij.vars.out,1,[]),vars_out(is_k)]});
        dom_ij = struct('out',[op_ij.dom.out;dom_out(is_k,:)]);
        if isa(op_ij,'tensopvar')
            % Nested operator: add the new variables with degree 1
            deps_ij = [op_ij.deps.out,ones(1,nnz(is_k))];
            Cop_out.ops{i,j} = lift_output(op_ij,vars_ij,dom_ij,[],deps_ij);
        else
            Cop_out.ops{i,j} = lift_output(op_ij,vars_ij,dom_ij);
        end
    end
end

% Update the variables, domains, and dependencies, on the output variables
% that are in the factors
if any(host(:))
    Cop_out.chkval = false;
    Cop_out = initialize(Cop_out);
end

% Add the variables that no factor maps to, in sorted order as for any
% tensopvar object
if numel(Cop_out.vars.out)~=M
    [vars_s,ord] = sort(vars_out);
    [~,idcs_s] = ismember(Cop_out.vars.out,vars_s);
    depmat_s = zeros(d,M);
    depmat_s(:,idcs_s) = Cop_out.depmat.out;
    Cop_out.vars.out = vars_s;
    Cop_out.dom.out = dom_out(ord,:);
    Cop_out.depmat.out = depmat_s;
end

% Lift a scalar output to a vector: repeating the scalar operator along the
% rows only repeats the first factor
if m_new~=size(Cop,1)
    Cop_out = repmat(Cop_out,m_new,1);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [vars_out,dom_out,dim_out] = parse_inputs(vars,dom,dims,vars0,dom0,matdim0)
% [VARS_OUT,DOM_OUT,DIM_OUT] = PARSE_INPUTS(VARS,DOM,DIMS,VARS0,DOM0,
% MATDIM0) checks the inputs specifying the output space to which the
% operator is lifted, and returns them in standard form.
%
% INPUTS
% - vars, dom: the output variables (cellstr) and their domains (M x 2
%           array) of the lifted operator, or structs with field 'out'
%           containing these;
% - dims:   [], or the number of rows of the lifted operator, specified as
%           a scalar, as 1 x 2 array [rows,columns], or as struct with
%           field 'out';
% - vars0, dom0: structs with field 'out', the output variables and
%           domains of the operator that is lifted;
% - matdim0: 1 x 2 array of the matrix dimensions of that operator;
%
% OUTPUTS
% - vars_out: 1 x M cellstr of output variables of the lifted operator;
% - dom_out: M x 2 'double' array of the associated domains;
% - dim_out: number of rows of the lifted operator: matdim0(1), or dims if
%           specified, in which case matdim0(1) must be 1 or equal dims;
%
% Check that the output variables are properly specified, as cellstr or as
% struct with field 'out'
if isa(vars,'struct') && isfield(vars,'out')
    vars = vars.out;
end
if ~iscellstr(vars)
    error("Output variable names should be specified as 'cellstr' object, or struct with field 'out'.")
end
vars_out = reshape(vars,1,[]);          M = numel(vars_out);
if numel(unique(vars_out))~=M
    error("Output variable names must be unique.")
end

% Check that the domains are properly specified, as numeric array or as
% struct with field 'out'
if isa(dom,'struct') && isfield(dom,'out')
    dom = dom.out;
end
if ~isnumeric(dom) || ~all(size(dom)==[M,2])
    error("Number of output domains should match number of output variables, and be specified as M x 2 array or struct with field 'out'.")
end
dom_out = dom;

% The output variables of the operator must be kept, on the same domain,
% as lifting extends the output function space but does not restrict it
vars_old = reshape(vars0.out,1,[]);
[is_old,idcs_old] = ismember(vars_old,vars_out);
if ~all(is_old)
    error("The output variables of the operator must be among the output variables to lift to.")
end
if ~isempty(vars_old) && any(any(dom_out(idcs_old,:)~=dom0.out & ~(isnan(dom_out(idcs_old,:)) & isnan(dom0.out))))
    error("The domains of the output variables of the operator cannot be changed.")
end

% Check that the number of rows is properly specified
dim_out = matdim0(1);
if ~isempty(dims)
    if isa(dims,'struct') && isfield(dims,'out')
        dims = dims.out;
    end
    if ~isnumeric(dims) || ~(isscalar(dims) || numel(dims)==2) || any(dims~=round(dims)) || any(dims<0)
        error("Dimensions must be specified as nonnegative integer, 1x2 array of nonnegative integers, or struct with field 'out'.")
    elseif numel(dims)==2 && dims(2)~=matdim0(2)
        error("Number of columns cannot be changed when lifting the output.")
    end
    dim_out = dims(1);
    % Only a scalar output can be lifted to a vector: all components of
    % the lifted output are then equal to this scalar
    if dim_out~=matdim0(1) && matdim0(1)~=1
        error("Only a scalar-valued output can be lifted to a vector-valued output.")
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
