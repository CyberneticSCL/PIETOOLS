function Cop = initialize(Cop)
% COP = INITIALIZE(COP) takes a 'tensopmat' object representing a block
% matrix of tensor-PI operators, and fills in the empty fields based on the
% information on the individual tensor-PI operators in each element.
%
% INPUTS
% - Cop:    'tensopmat' object with field Cop.ops specified as a p x d cell
%           array of 'tensopvar' or 'intop' objects, representing a
%           concatenation of tensor-PI operators, or of matrices of type
%           'double', 'polynomial' or 'dpvar' acting on R (see
%           `tensopmat`). The domains of the variables in such matrices are
%           taken from the other blocks, or from Cop.vars and Cop.dom if
%           already set;
%
% OUTPUTS
% - Cop:    'tensopmat' object representing the same tensor-PI operator as
%           the input, but with fields 'dims', 'vars', 'dom', and 'depmat'
%           filled in based on the values of Cop.ops;
%
% NOTES
% - Empty blocks Cop.ops{i,j} represent zero operators.
% - Blocks of type 'opvar', 'opvar2d' or 'sopvar' are converted to
%   'tensopvar' objects in Cop.ops.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - initialize
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
% DJ, 10/15/2026: Initial coding;

% Determine the input domain from the first row of operators
[p,d] = size(Cop.ops);
dims_in = zeros(d,1);   vars_in = cell(1,0);    dom_in = zeros(0,2);
is_nat = native_mask(Cop.ops);        % blocks of type double/polynomial/dpvar
for j=1:d
    for i=1:p
        op_ij = Cop.ops{i,j};
        if isempty(op_ij)
            if i==p
                error("Input dimensions of the operator are ambiguous.")
            end
            continue
        end
        % Get the column dimensions of the operators
        dims_in(j) = size(op_ij,2);
        if is_nat(i,j)
            break       % acts on R: no input variables
        end

        % Add the input and output variables to the unique list
        vars_j = op_ij.vars;                       
        vars_in_full = [vars_in,vars_j.in];
        [vars_in,idcs1_in] = unique(vars_in_full);
    
        % Add the associated domains to the unique list
        dom_j = op_ij.dom;
        dom_in_full = [dom_in; dom_j.in];
        dom_in = dom_in_full(idcs1_in,:);
        break
    end
end

% Determine the output domain from the first column of operators
dims_out = zeros(p,1);   vars_out = cell(1,0);    dom_out = zeros(0,2);
vars_rows = cell(p,1);      % output variables of each row
is_nat_row = any(is_nat,2) & ~any(~is_nat & ~cellfun('isempty',Cop.ops),2);
for i=1:p
    if is_nat_row(i)
        % Row of blocks of type double/polynomial/dpvar only: the output variables are those
        % appearing in the blocks, with domains from elsewhere
        j = find(is_nat(i,:),1);
        dims_out(i) = size(Cop.ops{i,j},1);
        vars_i = cell(1,0);
        for j=find(is_nat(i,:))
            vars_i = [vars_i,native_vars(Cop.ops{i,j})];
        end
        vars_i = unique(vars_i);
        vars_rows{i} = vars_i;
        dom_i = lookup_dom(Cop,is_nat,vars_i);
        vars_out_full = [vars_out,vars_i];
        [vars_out,idcs1_out] = unique(vars_out_full);
        dom_out_full = [dom_out; dom_i];
        dom_out = dom_out_full(idcs1_out,:);
        continue
    end
    for j=1:d
        op_ij = Cop.ops{i,j};
        if isempty(op_ij)
            if j==d
                error("Output dimensions of the operator are ambiguous.")
            end
            continue
        end
        % Get the row dimensions of the operators
        dims_out(i) = size(op_ij,1);
        if is_nat(i,j)
            continue    % the other blocks in the row set the output variables
        end

        % Add the input and output variables to the unique list
        vars_i = op_ij.vars;
        vars_rows{i} = vars_i.out;
        vars_out_full = [vars_out,vars_i.out];
        [vars_out,idcs1_out] = unique(vars_out_full);

        % Add the associated domains to the unique list
        dom_i = op_ij.dom;
        dom_out_full = [dom_out; dom_i.out];
        dom_out = dom_out_full(idcs1_out,:);
        break
    end
end

% Check how the input and output domains of the operators relate to those
% of the full block
N = numel(vars_in);     M = numel(vars_out);
depmat_in = zeros(d,N);
depmat_out = zeros(p,M);
for i=1:p
    for j=1:d
        % Check that the operator is of appropriate type
        op_ij = Cop.ops{i,j};
        if isempty(op_ij)
            continue
        end
        if is_nat(i,j)
            % Block of type double/polynomial/dpvar: acts on R, so column j has depmat.in=0
            if dims_in(j)~=size(op_ij,2)
                error("Column dimensions of operators in same column-block must match.")
            elseif dims_out(i)~=size(op_ij,1)
                error("Row dimensions of operators in same row-block must match.")
            elseif ~all(ismember(native_vars(op_ij),vars_rows{i}))
                error("Variables in the block of type 'double', 'polynomial' or 'dpvar' on row %d, column %d must be output variables of the other blocks in the row.",i,j)
            end
            if is_nat_row(i)
                depmat_out(i,ismember(vars_out,vars_rows{i})) = 1;
            end
            continue
        end
        if isa(op_ij,'opvar') || isa(op_ij,'opvar2d')
            op_ij = tensopvar(opvar2sopvar(op_ij));
            Cop.ops{i,j} = op_ij;
        elseif isa(op_ij,'sopvar')
            op_ij = tensopvar(op_ij);
            Cop.ops{i,j} = op_ij;
        elseif ~isa(op_ij,'tensopvar') && ~isa(op_ij,'intop')
            error("Blocks in 'tensopmat' structure must be of type 'tensopvar', 'intop', 'double', 'polynomial' or 'dpvar'.")
        end
        if isa(op_ij,'tensopvar') && ~op_ij.chkval
            op_ij = initialize(op_ij);
        elseif isa(op_ij,'intop') && ~op_ij.chkval && ~isvalid(op_ij)
            error(['Operator P.ops{',num2str(i),',',num2str(j),'} is not properly declared.'])
        end

        % Check that the dimensions match
        if dims_in(j)~=size(op_ij,2)
            error("Column dimensions of operators in same column-block must match.")
        end
        if dims_out(i)~=size(op_ij,1)
            error("Row dimensions of operators in same row-block must match.")
        end

        % Check that the input variables match
        vars_ij = op_ij.vars;      deps_ij = op_ij.deps;
        if any(is_nat(:,j)) && any(deps_ij.in)
            error("Blocks of type 'double', 'polynomial' or 'dpvar' act on R, so cannot share column %d with an operator acting on a function space.",j)
        end
        [~,idcs1_in,idcs2_in] = unique([vars_in,vars_ij.in],'stable');
        if numel(idcs1_in)>N
            error("Operators in the same column should be defined in terms of the same input variables.")
        end
        idcs2_in = idcs2_in(N+1:end);
        depmat_in(j,idcs2_in) = deps_ij.in;
        % Check that the output variables match
        [~,idcs1_out,idcs2_out] = unique([vars_out,vars_ij.out],'stable');
        if numel(idcs1_out)>M
            error("Operators in the same row should be defined in terms of the same output variables.")
        end
        idcs2_out = idcs2_out(M+1:end);
        depmat_out(i,idcs2_out) = deps_ij.out;

        % Check that the domains match
        dom_ij = op_ij.dom;
        if ~isequal(dom_in(idcs2_in,:),dom_ij.in)
            error("Identical input variables should have identical domains.")
        end
        if ~isequal(dom_out(idcs2_out,:),dom_ij.out)
            error("Identical output variables should have identical domains.")
        end
    end
end

% Set the input and output domains of the operators
dims = struct();        dims.in = dims_in;      dims.out = dims_out;
vars = struct();        vars.in = vars_in;      vars.out = vars_out;
dom = struct();         dom.in = dom_in;        dom.out = dom_out;
depmat = struct();      depmat.in = depmat_in;  depmat.out = depmat_out;

Cop.dims = dims;
Cop.vars = vars;
Cop.dom = dom;
Cop.depmat = depmat;
Cop.chkval = true;

end


function dom = lookup_dom(Cop,is_nat,names)
% DOM = LOOKUP_DOM(COP,IS_NAT,NAMES) returns the M x 2 domains of the variables NAMES (1 x M cellstr)
% that appear in the blocks of type 'double', 'polynomial' or 'dpvar' (marked in IS_NAT), which carry no
% domain themselves. These are taken from the other blocks, and from Cop.vars, Cop.dom if already set.
dom = zeros(numel(names),2);
if isempty(names)
    return
end
pool_vars = cell(1,0);      pool_dom = zeros(0,2);
for k=find(~is_nat(:)' & ~cellfun('isempty',Cop.ops(:))')
    op = Cop.ops{k};
    if isa(op,'tensopvar') || isa(op,'intop') || isa(op,'sopvar')
        vars_k = op.vars;       dom_k = op.dom;
        pool_vars = [pool_vars,vars_k.in,vars_k.out];
        pool_dom = [pool_dom; dom_k.in; dom_k.out];
    end
end
if numel(Cop.vars.in)==size(Cop.dom.in,1) && numel(Cop.vars.out)==size(Cop.dom.out,1)
    pool_vars = [pool_vars,Cop.vars.in,Cop.vars.out];
    pool_dom = [pool_dom; Cop.dom.in; Cop.dom.out];
end
[is_found,idcs] = ismember(names,pool_vars);
if ~all(is_found)
    error("The domain of variable '%s' in a block of type 'double', 'polynomial' or 'dpvar' cannot be determined: specify the variables and their domains as inputs 'vars' and 'dom' to the constructor, or in Cop.vars.out and Cop.dom.out.",names{find(~is_found,1)})
end
dom = pool_dom(idcs,:);
end
