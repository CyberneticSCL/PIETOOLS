function [S,Sd,Z3deg] = int_outer(Z1deg,Z2deg)
% Given degmats of Z1 = length(m) and Z2 = length(n)
% we compute int(Z1*Z2^T, t) = S(I_n\otimes Z3) = (I_m\otimes Z3^T)Sd.
% S is size m x (n*length(Z3))
% Sd is size (m*length(Z3)) x n.

Z1deg = Z1deg(:);  Z2deg = Z2deg(:); % make it column
m  = numel(Z1deg); n  = numel(Z2deg); % number of monomials in Z1,Z2

% Pairwise sums Z3_ij and mapping to unique degrees Z3deg
sumMat = Z1deg + Z2deg.';                    % m x n, degrees Z3_ij
[Z3deg, ~, idx] = unique(sumMat(:), 'sorted');
m3 = numel(Z3deg);
idxMat = reshape(idx, m, n);           % k(i,j), k in 1:length(Z3deg)

% Integrated basis degrees: Z3_ij + 1  -> Z3deg_int = Z3deg + 1
Z3deg_int = Z3deg + 1;

% Weight per unique product degree: w_k = 1/(Z3deg(k)+1)
w = 1./double(Z3deg_int);               % m3 x 1

% For each (i,j), the nonzero value is w_{k(i,j)}
vals = w(idxMat(:));                   % length m*n

% ----- Build S: m x (n*m3) -----
% Nonzeros at (row=i, col=(j-1)*m3 + k(i,j)) with value w(k(i,j))
rowsS = repmat((1:m).', n, 1);                         % length m*n
colShiftS = kron((0:n-1).'*m3, ones(m,1));             % length m*n
colsS = idxMat(:) + colShiftS;                         % length m*n
S = sparse(rowsS, colsS, vals, m, n*m3);

% ----- Build Sd: (m*m3) x n -----
% Nonzeros at (row=(i-1)*m3 + k(i,j), col=j) with value w(k(i,j))
[I,J] = ndgrid(1:m, 1:n);
rowsSd = (I(:)-1)*m3 + idxMat(:);                      % length m*n
colsSd = J(:);                                         % length m*n
Sd = sparse(rowsSd, colsSd, vals, m*m3, n);
end