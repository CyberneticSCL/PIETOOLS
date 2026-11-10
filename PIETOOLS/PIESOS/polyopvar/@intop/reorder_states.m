function Kop_out = reorder_states(Kop,state_order,dims,depmat)
% KOP_OUT = REORDER_STATES(KOP,STATE_ORDER,DIMS,DEPMAT) reorders the
% integrals defining the PI functional KOP, acting on a distributed
% monomial, to match the new order of the state variables in the
% distributed monomial, so that
%   KOP_OUT*(x_{i1} o ... o x_{id}) = KOP*(x1 o ... o xd)
% where o denotes the tensor product.
%
% INPUTS
% - Kop:            m x n 'intop' object representing a functional on a
%                   degree-d distributed monomial;
% - state_order:    1 x d integer array, with values ranging from 1 to d,
%                   specifying the new order of the d state variables in
%                   the degree-d distributed monomial;
% - dims:           d x 1 integer array specifying the vector size of each
%                   of the d state variables, such that
%                       size(Kop,2)==prod(dims);
% - depmat:         d x N logical array specifying for each of the N
%                   spatial variables in terms of which the functional is
%                   defined which state variables depend on this spatial
%                   variable, such that
%                       Kop.depmat.in = sum(depmat,1);
%
% OUTPUTS
% - Kop_out:        m x n 'intop' object representing the same functional
%                   as the input, acting on the reordered monomial;
%
% NOTES
% - If depmat has more than N columns, only the columns of the spatial
%   variables on which at least one state depends are retained.
% - This kind of reordering of the state vectors in the distributed
%   monomial is possible because they are vectors, so that
%   kron(x,y) = P*kron(y,x) for a suitable permutation matrix P. If x or y
%   is matrix-valued, another permutation matrix would be needed,
%   kron(x,y) = P1*kron(y,x)*P2.
% - The columns of the functional are permuted by permuting the blocks of
%   nZ columns (one for each of the n matrix columns) of the coefficient
%   matrices.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - reorder_states
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

% Check that the functional is properly specified
if nargin<=3
    error("Not enough inputs.")
end
if ~isa(Kop,'intop')
    error("Functional must be specified as 'intop' object.")
end
if ~Kop.chkval && ~isvalid(Kop)
    error("The functional is not properly declared.")
end
varnames = Kop.varnames;
N = numel(varnames.in_dum);
degs = Kop.depmat.in;

% Check that the variable order is properly specified
if ~isnumeric(state_order)
    error("Order of variables must be specified as 1 x d array of integers.")
elseif numel(unique(state_order))~=numel(state_order)
    error("A new order must be specified for each of the factors in the distributed monomial.")
end
state_order = reshape(state_order,1,[]);
d = numel(state_order);
if all(state_order==1:d)
    % No reordering is performed
    Kop_out = Kop;
    return
end

% Check that the dimensions are properly specified
if ~isnumeric(dims) || numel(dims)~=d
    error("A vector dimension must be specified for each of the factors in the distributed monomials.")
elseif prod(dims)~=size(Kop,2)
    error("Vector dimension of the distributed monomial does not match column dimension of the functional.")
end

% Check that the dependency array is appropriate
if size(depmat,2)>N
    % Retain only relevant columns of the dependency array
    depmat = depmat(:,any(depmat,1));
end
if ~all(size(depmat)==[d,N]) || ~islogical(depmat)
    error("Dependency array must be specified as d x N logical array, for d factors and N spatial variables.")
elseif any(sum(depmat,1)~=degs)
    error("Input degrees of variables in functional do not match specified dependency array.")
end

% Determine the new order of the spatial variables
var_order = cell(1,N);
for i=1:N
    var_idcs = zeros(1,d);
    var_idcs(depmat(:,i)) = 1:degs(i);      % index of variable i in each monomial factor
    var_idcs = var_idcs(state_order);
    deps_i = depmat(state_order,i);
    var_order{i} = var_idcs(deps_i);
end

% Reorder the spatial (dummy) variables in the functional
Kop = reorder_vars(Kop,var_order);

% Reorder the columns of the parameters to match new order of factors
%   kron(x1,x2) = Pmat*kron(x2,x1)
if size(Kop,2)==1
    % No reordering in scalar case
    Kop_out = Kop;
    return
end
% For each column of the functional, establish which element of each of the
% factors in the monomial it corresponds to
mon_idx_old = zeros(1,0);
for i=1:d
    mon_idx_old = [repelem(mon_idx_old,dims(i),1),repmat((1:dims(i))',size(mon_idx_old,1),1)];
end
% Switch the order of the factors in the monomial
mon_idx_tmp = mon_idx_old(:,state_order);   % column i in Kop is still linked with row i in mon_idx_tmp
% Swap Kop-->Kop(:,idcs) so that mon_idx_tmp(idcs,:) matches expected order
[~,idcs] = sortrows_integerTable(mon_idx_tmp);
% Each matrix column corresponds to a block of nZ columns of the coefficients
[~,~,nZ] = Zdims(Kop.Z.in);
colidx = reshape((idcs(:)'-1)*nZ + (1:nZ)',1,[]);
params = Kop.params;
for k=1:numel(params)
    params{k} = params{k}(:,colidx);
end
Kop_out = Kop;
Kop_out.params = params;
Kop_out.chkval = true;

end
