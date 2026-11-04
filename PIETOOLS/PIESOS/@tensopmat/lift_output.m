function Cop_out = lift_output(Cop,vars,dom,dims,deps)
% COP_OUT = LIFT_OUTPUT(COP,VARS,DOM,DIMS,DEPS) returns the 'tensopmat'
% object COP_OUT representing the block operator COP with the output of
% each row-block lifted to the function space defined by VARS, DOM and
% DEPS: the output is extended with variables it was not yet defined in,
% along which it is constant. See also 'lift_output' for 'tensopvar'
% objects.
%
% INPUTS
% - Cop:    p x d 'tensopmat' object representing an operator
%               Cop: L2^n[dom_in] --> L2^m[dom_out];
% - vars:   1 x M cellstr specifying the names of the output variables of
%           the lifted operator. Must contain all output variables of Cop.
%           Can also be specified as struct with field 'out' containing
%           this cellstr (other fields, e.g. 'in', are ignored);
% - dom:    M x 2 'double' array specifying the intervals on which these
%           variables are defined. Must match the domains of the output
%           variables of Cop. Can also be specified as struct with field
%           'out' (other fields are ignored);
% - dims:   (optional) p x 1 array of integers (or struct with field
%           'out') specifying the row dimension of each row-block of the
%           lifted operator, or a scalar for the same dimension of all
%           row-blocks. Defaults to the dimensions Cop.dims.out. A
%           scalar-valued row-block can be lifted to a vector-valued one,
%           repeating the scalar output in each row. If Cop has a single
%           (scalar-valued) row-block, dims can be a q x 1 array, to repeat
%           the row-block q times, with dimension dims(i) for the ith copy;
% - deps:   (optional) p x M array of nonnegative integers (or struct with
%           field 'out') specifying the degree of each row-block in each of
%           the output variables vars, i.e. the number of copies of the
%           variable in the output function space of the row-block (see
%           Cop.depmat.out). A 1 x M array specifies the same degrees for
%           all row-blocks. Must be at least the degrees of Cop, and
%           defaults to these degrees, or to 1 for variables that a
%           row-block is not yet defined in (so that all row-blocks are
%           lifted to the same function space);
%
% OUTPUTS
% - Cop_out: p x d 'tensopmat' object representing the lifted operator,
%           with
%               Cop_out.vars.out = sort(vars);
%               Cop_out.depmat.out = deps (sorted accordingly);
%
% NOTES
% - Blocks of type 'tensopvar' are lifted by 'lift_output', which
%   determines in which factors the new copies of the variables are
%   introduced.
% - Blocks of type 'double', 'polynomial' or 'dpvar' acting on R (see
%   'tensopmat') are not changed (a constant is a constant function),
%   apart from repeating the rows to lift a scalar to a vector. This also
%   holds for a row-block consisting of such blocks only, which maps to
%   the variables that appear in the blocks, and for which the output
%   space is just set to the variables vars: this is possible for degree
%   at most 1 in each variable.
% - Blocks of type 'intop' (functionals) cannot be lifted, so a row-block
%   containing such a block must not be lifted along any variable.
% - A new variable with degree 0 in each row-block is only added to the
%   list of variables, without changing the operator.
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

% Check that the block operator is properly declared
if ~isa(Cop,'tensopmat')
    error("First argument must be of type 'tensopmat'.")
end
if ~Cop.chkval
    Cop = initialize(Cop);
end
[p,d] = size(Cop.ops);
if p==0 || d==0
    error("Cannot lift an operator without blocks.")
end
if nargin<=3
    dims = [];
end
[vars_out,dom_out,dims_out] = parse_inputs(vars,dom,dims,Cop.vars,Cop.dom,Cop.dims.out);
M = numel(vars_out);

% Lifting a scalar to a vector with several row-blocks: repeat the single
% row-block, to lift each copy to the dimension of the row-block
if numel(dims_out)>p
    Cop.ops = repmat(Cop.ops,numel(dims_out),1);
    Cop.dims.out = repmat(Cop.dims.out,numel(dims_out),1);
    Cop.depmat.out = repmat(Cop.depmat.out,numel(dims_out),1);
    p = numel(dims_out);
end

% Extract the degrees of each row-block in the new output variables, with
% degree 0 for the variables that Cop is not defined in
[~,idcs_old] = ismember(Cop.vars.out,vars_out);
deps_old = zeros(p,M);              deps_old(:,idcs_old) = Cop.depmat.out;

% Check that the degrees of the lifted operator are properly specified
if nargin<=4 || isempty(deps)
    deps = max(deps_old,1);
elseif isa(deps,'struct') && isfield(deps,'out')
    deps = deps.out;
end
if ~isnumeric(deps) || any(deps(:)~=round(deps(:))) || any(deps(:)<0)
    error("Degrees must be specified as array of nonnegative integers.")
elseif isequal(size(deps),[1,M])
    deps = repmat(deps,p,1);
elseif ~isequal(size(deps),[p,M])
    error("Degrees must be specified as 1 x M or p x M array, with p the number of row-blocks and M the number of output variables.")
end
if any(any(deps<deps_old))
    error("The degrees of the lifted operator cannot be smaller than those of the operator.")
end

% Lift the blocks in each row-block that changes
is_nat = native_mask(Cop.ops);
Cop_out = Cop;
for i=1:p
    is_changed = any(deps(i,:)~=deps_old(i,:));
    if ~is_changed && dims_out(i)==Cop.dims.out(i)
        continue
    end
    deps_i = deps(i,:);
    in_i = deps_i>0;
    nrep = dims_out(i)/max(Cop.dims.out(i),1);          % repetitions, to lift scalar to vector
    if any(is_nat(i,:)) && all(is_nat(i,:) | cellfun('isempty',Cop.ops(i,:))) && any(deps_i>1)
        error("A row-block consisting of blocks of type 'double', 'polynomial' or 'dpvar' can only be lifted to degree at most 1 in each variable.")
    end
    for j=1:d
        op_ij = Cop.ops{i,j};
        if isempty(op_ij)
            continue
        elseif is_nat(i,j)
            % Matrix acting on R: a constant is a constant function, so
            % only repeat the matrix to lift a scalar to a vector
            if nrep~=1
                Cop_out.ops{i,j} = repmat(op_ij,nrep,1);
            end
        elseif isa(op_ij,'tensopvar')
            % Lift to the variables on which the row-block is defined
            Cop_out.ops{i,j} = lift_output(op_ij,vars_out(in_i),dom_out(in_i,:),dims_out(i),deps_i(in_i));
        else
            % Functional, with no output variables
            if is_changed
                error("Blocks of type 'intop' (functionals) cannot be lifted along output variables (row-block %d, column-block %d).",i,j)
            end
            % Repeat the rows of the functional (vertcat returns a tensopmat)
            if nrep~=1
                Cop_out.ops{i,j} = ones(nrep,1)*op_ij;
            end
        end
    end
end

% Update the dimensions, variables, domains and degrees
Cop_out.chkval = false;
Cop_out = initialize(Cop_out);
% Set the output space to that specified. For the row-blocks consisting
% only of matrices acting on R, and for variables that no block maps to,
% this is not determined by the blocks
[vars_s,ord] = sort(vars_out);
Cop_out.vars.out = vars_s;
Cop_out.dom.out = dom_out(ord,:);
Cop_out.depmat.out = deps(:,ord);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [vars_out,dom_out,dims_out] = parse_inputs(vars,dom,dims,vars0,dom0,dims0)
% [VARS_OUT,DOM_OUT,DIMS_OUT] = PARSE_INPUTS(VARS,DOM,DIMS,VARS0,DOM0,
% DIMS0) checks the inputs specifying the output space to which the
% operator is lifted, and returns them in standard form.
%
% INPUTS
% - vars, dom: the output variables (cellstr) and their domains (M x 2
%           array) of the lifted operator, or structs with field 'out'
%           containing these;
% - dims:   [], or the row dimensions of the row-blocks of the lifted
%           operator, as scalar, p x 1 array (or q x 1 array if p=1), or
%           struct with field 'out';
% - vars0, dom0: structs with field 'out', the output variables and
%           domains of the operator that is lifted;
% - dims0:  p x 1 array of the row dimensions of its row-blocks;
%
% OUTPUTS
% - vars_out: 1 x M cellstr of output variables of the lifted operator;
% - dom_out: M x 2 'double' array of the associated domains;
% - dims_out: array of the row dimensions of the row-blocks of the lifted
%           operator (q x 1, with q>p only if p=1, to repeat the single
%           row-block): dims0, or dims if specified, in which case each
%           element must equal that of dims0 or the latter must be 1;
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

% Check that the row dimensions are properly specified
p = numel(dims0);
dims_out = reshape(dims0,[],1);
if ~isempty(dims)
    if isa(dims,'struct') && isfield(dims,'out')
        dims = dims.out;
    end
    if ~isnumeric(dims) || ~(isscalar(dims) || numel(dims)==p || p==1) || any(dims(:)~=round(dims(:))) || any(dims(:)<0)
        error("Dimensions must be specified as scalar or array of nonnegative integers, with an element for each row-block, or as struct with field 'out'.")
    end
    dims_out = reshape(dims,[],1).*ones(p,1);
    % Only a scalar output can be lifted to a vector: all components of
    % the lifted output are then equal to this scalar
    if any(dims_out~=dims0(:) & dims0(:)~=1) || (numel(dims_out)>p && any(dims0(:)~=1))
        error("Only a scalar-valued output can be lifted to a vector-valued output.")
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
