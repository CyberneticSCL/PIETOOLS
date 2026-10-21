function [Z,M] = merge_monomial_product(Z1,Z2)
% [Z,M] = MERGE_MONOMIAL_PRODUCT(Z1,Z2) collapses the elementwise product of
% two monomial bases in the SAME set of variables into a single monomial
% basis, returning the map between them. That is, it returns Z and M such
% that
%
%       kron(Z1(x),Z2(x)) = M*Z(x)
%
% where, following the sopvar convention,
%
%       Zi(x) = kron(Zi{1}(x_1),...,Zi{N}(x_N)).
%
% This is needed whenever a coefficient matrix carries monomials in x on
% both sides of an operation that generates further monomials in x, as
% happens when the integration variable of a positive operator is
% eliminated: the pre-existing basis and the generated basis must be merged
% before the result can be stored as a single sopvar/sdopvar parameter.
%
% INPUTS
% - Z1, Z2: 1 x N cell arrays of integer column vectors, representing
%           monomial bases in the same N variables
%
% OUTPUTS
% - Z:      1 x N cell array of integer column vectors, with Z{i} the set of
%           degrees obtainable as a sum of a degree in Z1{i} and a degree in
%           Z2{i}
% - M:      (n1*n2) x n sparse 0/1 matrix with exactly one nonzero per row,
%           where ni = prod(cellfun(@numel,Zi)) and n = prod(cellfun(@numel,Z))
%
% MP, 08/22/2026: Initial coding

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - merge_monomial_product
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

Z1 = Z1(:).';
Z2 = Z2(:).';
N = numel(Z1);
if numel(Z2)~=N
    error("Monomial bases must be defined in the same number of variables.")
end

% Degenerate case: no variables, so both bases are the scalar 1.
if N==0
    Z = cell(1,0);
    M = sparse(1,1,1,1,1);
    return
end

% Degrees appearing in the product, per variable.
Z = cell(1,N);
for i=1:N
    d1 = reshape(Z1{i},[],1);
    d2 = reshape(Z2{i},[],1);
    Z{i} = unique(reshape(d1+d2.',[],1));
end

% Build the full degree matrices, using the convention that variable 1 is
% the outermost (slowest) index of the Kronecker product.
degmat_1 = full_degmat(Z1);
degmat_2 = full_degmat(Z2);
degmat   = full_degmat(Z);

n1 = size(degmat_1,1);
n2 = size(degmat_2,1);
n  = size(degmat,1);

% Row (u-1)*n2+v of kron(Z1(x),Z2(x)) is the monomial of degree
% degmat_1(u,:)+degmat_2(v,:).
degmat_prod = repelem(degmat_1,n2,1) + repmat(degmat_2,n1,1);

[tf,idx] = ismember(degmat_prod,degmat,'rows');
if ~all(tf)
    error("Internal error: product monomial not found in merged basis.")
end

M = sparse(1:(n1*n2),idx,ones(n1*n2,1),n1*n2,n);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function degmat = full_degmat(Zc)
% Expand a 1 x N cell of per-variable degree vectors into the degree matrix
% of the full tensor basis kron(Zc{1},...,Zc{N}), with variable 1 outermost.

degmat = zeros(1,0);
for i=1:numel(Zc)
    Zi = reshape(Zc{i},[],1);
    degmat = [kron(degmat,ones(numel(Zi),1)), kron(ones(size(degmat,1),1),Zi)];
end

end
