function F = monomial_shift_right(C, nz)
% Given C of size (m*nz)-by-n, computes F such that
%
%   kron(I_m, Z') * C = F * kron(I_n, Z)
%
% C is (m*nz)-by-n
% F is m-by-(n*nz)
%
% This is the inverse of monomial_shift_sparse.

    [q, n] = size(C);

    if mod(q, nz) ~= 0
        error('Number of rows of C must be divisible by nz.');
    end

    m = q / nz;

    % Nonzero entries of shifted matrix C
    [row, j, val] = find(C);

    % Decode row:
    % row = (i-1)*nz + k
    i = floor((row - 1) / nz) + 1;     % original row index, 1,...,m
    k = mod(row - 1, nz) + 1;          % monomial index, 1,...,nz

    % Original column in F:
    % col = (j-1)*nz + k
    new_row = i;
    new_col = (j - 1) * nz + k;

    % Build sparse matrix in one shot
    F = sparse(new_row, new_col, val, m, n * nz);
end


% This should do the following MATH operations:
% C*kron(I_n, Z) = [C1*Z, C2*Z,..., Cn*Z]
% kron(I_m, Z')*FC = [Z'*FC1, ..., Z'*FCn]
% where FCi = vec(Ci')