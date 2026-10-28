function [Z,M] = merge_monomial_product(Z1,Z2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Z,M] = MERGE_MONOMIAL_PRODUCT(Z1,Z2) returns, for two monomial bases Z1,
% Z2 in the same N variables x, the monomial basis Z of all products of a
% monomial of Z1 and a monomial of Z2, and the map M between them:
%
%       kron(Z1(x),Z2(x)) = M*Z(x),
%
% where Zi(x) = kron(x_1.^Zi{1},...,x_N.^Zi{N}).
%
% INPUT
% - Z1, Z2: 'cell', 1 x N, of double column vectors of nonnegative integer
%           exponents: Zi{k} the degrees of x_k in basis Zi.
%
% OUTPUT
% - Z:      'cell', 1 x N, of double column vectors: Z{k} the degrees
%           obtainable as a sum of a degree in Z1{k} and a degree in Z2{k},
%           each once, in ascending order;
% - M:      'double', sparse, (n1*n2) x n, 0/1 with exactly one nonzero per
%           row, ni = prod(cellfun(@numel,Zi)) and
%           n = prod(cellfun(@numel,Z)); for N = 0, Z = cell(1,0) and M = 1.
%
% NOTES
% Z1 and Z2 with different numbers of entries are an error.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

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
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% MP, 08/22/2026: Initial coding
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), step comments in the body and INPUT and
%                OUTPUT lines for the local function; code unchanged. The
%                log entry above moved here from the help text, below the
%                licence block, which gains its last paragraph. The deleted
%                header text was the reason the merge is needed (the
%                elimination of an integration variable) and a reference
%                to the sopvar convention.

% Both bases as rows of N per-variable exponent vectors
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

% Numbers of monomials of Z1, Z2 and Z
n1 = size(degmat_1,1);
n2 = size(degmat_2,1);
n  = size(degmat,1);

% Row (u-1)*n2+v of kron(Z1(x),Z2(x)) is the monomial of degree
% degmat_1(u,:)+degmat_2(v,:).
degmat_prod = repelem(degmat_1,n2,1) + repmat(degmat_2,n1,1);

% Position of each product monomial in Z; M has its 1 there
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
% INPUT  Zc:     1 x N cell of double column vectors of exponents.
% OUTPUT degmat: prod(cellfun(@numel,Zc)) x N double, row t the exponents
%                of monomial t (1 x 0 for N = 0).

% Each variable appended as the innermost (fastest) index
degmat = zeros(1,0);
for i=1:numel(Zc)
    Zi = reshape(Zc{i},[],1);
    degmat = [kron(degmat,ones(numel(Zi),1)), kron(ones(size(degmat,1),1),Zi)];
end

end
