function Pmat = commat(m,n,opt)
% PMAT = COMMAT(M,N) returns the (M,N) commutation matrix PMAT, defined as
%   Pmat = [I_{m} o e_{1}^T;
%           :
%          [I_{m} o e_{n}^T]
% wher I_{m} is the m x m identity matrix, e_{i} the ith standard Euclidean
% basis vector in R^{n}, and o denotes the Kronecker product.
% The commutation matrix is so called because for A in R^{m x n} and B in
% R^{p x q}, we have
%       B o A = P_{m,p} (A o B) P_{n,q}^T
% Set opt='transpose' to return Pmat^T

r_idcs = (1:m)' + (0:n-1)*m;
c_idcs = (0:m-1)'*n + (1:n);
if nargin==3 && strcmp(opt,'transpose')
    Pmat = sparse(c_idcs,r_idcs,1,m*n,m*n);
else
    Pmat = sparse(r_idcs,c_idcs,1,m*n,m*n);
end

end

