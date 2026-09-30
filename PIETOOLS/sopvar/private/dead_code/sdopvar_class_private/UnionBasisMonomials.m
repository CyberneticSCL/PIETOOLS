function [Z, C1, C2] = UnionBasisMonomials(Z1, Z2)
% [Z,C1,C2] = UNIONBASISMONOMIALS(Z1,Z2) combines the monomial bases Z1 and
% Z2 into a unique basis Z, and return transformation matrices C1 and C2
% such that Z1 = C1*Z and Z2 = C2*Z
%
% INPUTS
% - Z1, Z2: 1 x N cell arrays of integer column vectors, representing
%           monomial vectors in N variables
%
% OUTPUTS
% - Z:      1 x N cell array of integer column vectors, corresponding to
%           the union of the column vectors in Z1 and Z2
% - C1,C2:  ni x n_new sparse matrix such that Zi(x) = Ci*Z(x), where Zi(x)
%           is the monomial vector defined by Zi, i.e.
%               Zi(x) = kron(Zi{1}(x),...,Zi{N}(x))
%
% MMP, SS, AT, DJ, 06/08/2026: Initial Coding

N = numel(Z1);
Z = cell(1,N);

% construct union basis
for i = 1:N
    Zi = union(Z1{i}, Z2{i});
    Z{i} = reshape(Zi,[],1);
end

% construct full monomial basis
degmat_1 = zeros(1,0);
degmat_2 = zeros(1,0);
degmat_new = zeros(1,0);
for dim = 1:N
    degmat_new = [kron(degmat_new, ones(size(Z{dim}))),  kron(ones(size(degmat_new,1), 1), Z{dim})];
    degmat_1 = [kron(degmat_1, ones(size(Z1{dim}))),  kron(ones(size(degmat_1,1), 1), Z1{dim})];
    degmat_2 = [kron(degmat_2, ones(size(Z2{dim}))),  kron(ones(size(degmat_2,1), 1), Z2{dim})];
end

% find identical rows 
[~, C1_index] =ismember(degmat_1,  degmat_new, 'rows');
[~, C2_index] =ismember(degmat_2,  degmat_new, 'rows');

n1 = size(degmat_1,1);
n2 = size(degmat_2,1);
n_new = size(degmat_new,1);
C1 = sparse(1:n1, C1_index, ones(n1, 1), n1, n_new);
C2 = sparse(1:n2, C2_index, ones(n2, 1), n2, n_new);

end