function Cop = initialize(Cop)
% COP = INITIALIZE(COP) takes a 'tensopvar' object representing a tensor-PI
% operator, and fills in the empty fields based on information on the PI
% operators defining each of the factors.
%
% INPUTS
% - Cop:    'tensopvar' object with field Cop.ops specified as a K x d cell
%           array of 'sopvar' (or 'intop' or 'tensopvar') objects,
%           representing a sum of K products of d PI operators;
%
% OUTPUTS
% - Cop:    'tensopvar' object representing the same tensor-PI operator as
%           the input, but with fields 'dims', 'vars', 'dom', and 'depmat'
%           filled in based on the values of Cop.ops;
%
% NOTES
% - Elements of Cop.ops of type 'opvar' or 'opvar2d' are converted to
%   'sopvar' objects.
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


% Extract the dimensions and domains of the factors defining the tensor-PI
% operator
[K,d] = size(Cop.ops);
dims = zeros(d,2);
vars = struct();        
dom = struct();
for i=1:K
    dims_in_i = zeros(d,1);     dims_out_i = zeros(d,1);
    vars_in_i = cell(1,0);      vars_out_i = cell(1,0);
    dom_in_i = zeros(0,2);      dom_out_i = zeros(0,2);
    for j=1:d
        % Check that the operator is of appropriate type
        op_ij = Cop.ops{i,j};
        if isa(op_ij,'opvar') || isa(op_ij,'opvar2d')
            op_ij = opvar2sopvar(op_ij);
            Cop.ops{i,j} = op_ij;
        elseif ~isa(op_ij,'sopvar') && ~isa(op_ij,'tensopvar') && ~isa(op_ij,'intop')
            error("Factors in tensor-PI operator must be of type 'sopvar'.")
        end

        % Get the dimensions of the operator
        [dim_out,dim_in] = size(op_ij);
        dims_out_i(j) = dim_out;
        dims_in_i(j) = dim_in;

        % Add the input and output variables to the unique list
        vars_j = op_ij.vars;                       
        vars_in_full = [vars_in_i,vars_j.in];  
        vars_out_full = [vars_out_i,vars_j.out];   
        [vars_in_i,idcs1_in,idcs2_in] = unique(vars_in_full);
        [vars_out_i,idcs1_out,idcs2_out] = unique(vars_out_full);

        % Add the associated domains to the unique list
        dom_j = op_ij.dom;
        dom_in_full = [dom_in_i; dom_j.in];     
        dom_out_full = [dom_out_i; dom_j.out];
        dom_in_i = dom_in_full(idcs1_in,:);
        dom_out_i = dom_out_full(idcs1_out,:);

        % Make sure the same variables have the same domains
        if any(any(dom_in_i(idcs2_in,:)~=dom_in_full)) ||...
            any(any(dom_out_i(idcs2_out,:)~=dom_out_full))
            error("Identical variables in different factors must have the same domains.")
        end
    end
    if i==1
        % Set the dimensions and domains
        dims(:,2) = dims_in_i;      dims(:,1) = dims_out_i;
        vars.in = vars_in_i;        vars.out = vars_out_i;
        dom.in = dom_in_i;          dom.out = dom_out_i;
    else
        if any(dims(:,2)~=dims_in_i) || any(dims(:,1)~=dims_out_i)
            error("Matrix dimensions of factors in each term must match.")
        elseif ~isequal(vars_in_i,vars.in) || ~isequal(vars_out_i,vars.out)
            error("Input and output variables of each term in the tensor-PI operator must match.")
        elseif ~isequal(dom_in_i,dom.in) || ~isequal(dom_out_i,dom.out)
            error("Domains of input and output variables of each term in the tensor-PI operator must match.")
        end
    end
end

% Check that the matrix dimensions make sense for elementwise product
type = Cop.otype;
if type(1)==1 && min(dims(:,1))~=max(dims(:,1))
    error("Row dimensions of factors must match for elementwise product.")
elseif type(2)==1 && min(dims(:,2))~=max(dims(:,2))
    error("Column dimensions of factors must match for elementwise product.")
end


% Set the dependency arrays
depmat = struct();
N = numel(vars.in);     M = numel(vars.out);
for i=1:K
    depmat_in_i = zeros(d,N);       depmat_out_i = zeros(d,M);
    for j=1:d
        % Check to what degree the operator depends on its variables
        op_ij = Cop.ops{i,j};
        vars_j = op_ij.vars;
        if isa(op_ij,'tensopvar')
            deps_j = op_ij.deps;
            deps_in_j = deps_j.in;
            deps_out_j = deps_j.out;
        else
            deps_in_j = ones(1,numel(vars_j.in));
            deps_out_j = ones(1,numel(vars_j.out));
        end

        % Set dependency in full list of variables
        [~,ia_in,ib_in] = intersect(vars.in,vars_j.in,'stable');
        depmat_in_i(j,ia_in) = deps_in_j(ib_in);
        [~,ia_out,ib_out] = intersect(vars.out,vars_j.out,'stable');
        depmat_out_i(j,ia_out) = deps_out_j(ib_out);
    end
    if i==1
        % Set the dependency array
        depmat.in = depmat_in_i;    depmat.out = depmat_out_i;
    else
        if ~isequal(depmat_in_i,depmat.in) || ~isequal(depmat_out_i,depmat.out)
            error("The same factor in different terms in the tensor-PI operator must depend on the same variables.")
        end
    end
end

% Declare the output object
Cop.vars = vars;
Cop.dom = dom;
Cop.dims = dims;
Cop.depmat = depmat;
Cop.chkval = true;

end