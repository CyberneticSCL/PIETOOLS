function FC = monomial_shift_left(C, nz)
% Given C of size m-by-(n*nz), computes FC such that
%
%   C * kron(I_n, Z) = kron(I_m, Z') * FC
%
% C  is m-by-(n*nz)
% FC is (m*nz)-by-n
%
% This avoids repeated sparse indexing assignments.

    [m, p] = size(C);

    if mod(p, nz) ~= 0
        error('Number of columns of C must be divisible by nz.');
    end

    n = p / nz;

    % Nonzero entries of C
    [i, col, val] = find(C);

    % Decode original column:
    % col = (j-1)*nz + k
    j = floor((col - 1) / nz) + 1;     % block index, 1,...,n
    k = mod(col - 1, nz) + 1;          % monomial index, 1,...,nz

    % New row in FC corresponds to pair (i,k):
    % row = (i-1)*nz + k
    new_row = (i - 1) * nz + k;
    new_col = j;

    % Build sparse matrix in one shot
    FC = sparse(new_row, new_col, val, m * nz, n);
end

% This should do the following MATH operations:
% C*kron(I_n, Z) = [C1*Z, C2*Z,..., Cn*Z]
% kron(I_m, Z')*FC = [Z'*FC1, ..., Z'*FCn]
% where FCi = vec(Ci')