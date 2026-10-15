function Pop = mat2tensopmat(Pmat,dims,vars,dom,depmat)
% POP = MAT2TENSOPMAT(PMAT,DIMS,VARS,DOM,DEPMAT) takes a matrix PMAT and
% returns a 'tensopmat' object POP representing the associated multiplier
% operator between the function spaces defined by DOM.
%
% INPUTS
% - Pmat:   m x n 'double' array, representing a multiplier operator;
% - dims:   (optional) struct with fields
%     in:   q x 1 array of integers adding up to n, specifying the vector
%           sizes of the inputs to the operator;
%    out:   p x 1 array of integers adding up to m, specifying the vector
%           sizes of the outputs to the operator. Defaults to a single
%           input and output block, dims.in = n and dims.out = m;
% - vars:   (optional) struct with fields
%     in:   1 x N cellstr array specifying the names of the input
%           variables to the operator;
%    out:   1 x M cellstr array specifying the names of the output
%           variables to the operator. Defaults to no variables;
% - dom:    (optional) struct with fields
%    in:    N x 2 'double' array specifying the intervals on which the
%           input variables are defined;
%   out:    M x 2 'double' array specifying the intervals on which the
%           output variables are defined. Defaults to no variables;
% - depmat: (optional) struct with fields
%    in:    q x N logical array specifying for each of the q inputs on
%           which of the N input variables it depends;
%   out:    p x M logical array specifying for each of the p outputs on
%           which of the M output variables it depends. Defaults to
%           no dependence on any variable;
%
% OUTPUTS
% - Pop:    'tensopmat' object with p x q blocks of type 'tensopvar',
%           representing the multiplier operator defined by Pmat;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mat2tensopmat
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

% Check that the matrix is properly specified
if ~isa(Pmat,'double')
    error("Input must be of type 'double'.")
end

% Check that the dimensions are appropriate
if nargin<=1
    dims = struct();
    dims.in = size(Pmat,2);     dims.out = size(Pmat,1);
end
if ~isa(dims,'struct') || ~isfield(dims,'in') || ~isfield(dims,'out')
    error("Dimensions must be specified as struct with fields 'in' and 'out'.")
end
if ~isnumeric(dims.in) || ~isnumeric(dims.out)
    error("Input and output vector dimensions should be specified as q x 1 integer arrays.")
end
dims.in = dims.in(:);       dims.out = dims.out(:);
q = numel(dims.in);         p = numel(dims.out);
if sum(dims.in)~=size(Pmat,2) || sum(dims.out)~=size(Pmat,1)
    error("Input and output dimensions must match the column and row dimensions of the matrix.")
end
nnr = [0;cumsum(dims.out)];     nnc = [0;cumsum(dims.in)];

% Check that the variables are properly specified
if nargin<=2
    vars = struct();
    vars.in = cell(1,0);        vars.out = cell(1,0);
end
if ~isa(vars,'struct') || ~isfield(vars,'in') || ~isfield(vars,'out')
    error("Variables must be specified as struct with fields 'in' and 'out'.")
end
if ~iscellstr(vars.in) || ~iscellstr(vars.out)
    error("Input and output variable names should be specified as 'cellstr' objects.")
end
vars.in = reshape(vars.in,1,[]);        N = numel(vars.in);
vars.out = reshape(vars.out,1,[]);      M = numel(vars.out);

% Check that the domains are properly specified
if nargin<=3
    dom = struct();
    dom.in = zeros(0,2);        dom.out = zeros(0,2);
end
if ~isa(dom,'struct') || ~isfield(dom,'in') || ~isfield(dom,'out')
    error("Domains must be specified as struct with fields 'in' and 'out'.")
end
if ~isnumeric(dom.in) || ~all(size(dom.in)==[N,2])
    error("Number of input domains should match number of input variables.")
end
if ~isnumeric(dom.out) || ~all(size(dom.out)==[M,2])
    error("Number of output domains should match number of output variables.")
end

% Check that the dependency arrays are appropriate
if nargin<=4
    depmat = struct();
    depmat.in = zeros(q,N);        depmat.out = zeros(p,M);
end
if ~isa(depmat,'struct') || ~isfield(depmat,'in') || ~isfield(depmat,'out')
    error("Dependency information must be specified as struct with fields 'in' and 'out'.")
end
if (~isnumeric(depmat.in) && ~islogical(depmat.in)) || ~all(size(depmat.in)==[q,N])
    error("Input dependency array must be of size q x N, for q inputs and N input variables.")
end
if (~isnumeric(depmat.out) && ~islogical(depmat.out)) || ~all(size(depmat.out)==[p,M])
    error("Output dependency array must be of size p x M, for p output and M output variables.")
end
if any(any(depmat.in>1)) || any(any(depmat.out>1))
    error("Conversion to higher-degree tensor-PI operators is currently not supported.")
end

% Build the actual operator
op_cell = cell(p,q);
for i=1:p
    % Extract the rows of the matrix for this output
    Pmat_i = Pmat(nnr(i)+1:nnr(i+1),:);
    % Set variables for this output
    vars_ij = vars;
    dom_ij = dom;
    vars_ij.out = vars.out(logical(depmat.out(i,:)));
    dom_ij.out = dom.out(logical(depmat.out(i,:)),:);
    for j=1:q
        % Extract the matrix associated with this block
        Pmat_ij = Pmat_i(:,nnc(j)+1:nnc(j+1));
        % Set the associated input variables
        vars_ij.in = vars.in(logical(depmat.in(j,:)));
        dom_ij.in = dom.in(logical(depmat.in(j,:)),:);
        % Declare the operator
        op_cell{i,j} = tensopvar( mat2sopvar(Pmat_ij,vars_ij,dom_ij));
    end
end
Pop = tensopmat();
Pop.ops = op_cell;
Pop.dims = dims;
Pop.vars = vars;
Pop.dom = dom;
Pop.depmat = depmat;
Pop.chkval = true;

end