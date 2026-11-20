function C = times(A,B)
% C = TIMES(A,B) returns a 'quadPoly' object representing the elementwise
% product of A and B, i.e. C = A.*B;
%
% INPUTS
% - A:  m x n 'double' or 'quadPoly' object;
% - B:  p x q 'double' or 'quadPoly' object, with p=m or p=1 or m=1, and
%       q=n or q=1 or n=1;
%
% OUTPUTS
% - C:  max(m,p) x max(n,q) 'quadPoly' object representing the elementwise
%       product of A and B, with singleton dimensions expanded as in
%       MATLAB's implicit expansion;
%
% NOTES
% - For two 'quadPoly' objects, the product is computed as
%   C(:) = diag(A(:))*B(:), so the same restrictions as in 'mtimes' apply:
%   left variables in A cannot appear as right variables in B, or vice
%   versa.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - times
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
% DJ, 01/21/2026: Initial coding;


if ~(isnumeric(A) || isa(A,'quadPoly')) || ~(isnumeric(B) || isa(B,'quadPoly'))
    error("Elementwise multiplication is only supported for 'quadPoly' and numeric objects.")
end

% Check that the dimensions are compatible, and expand singleton
% dimensions to the size of the product
[m,n] = size(A);    [p,q] = size(B);
if (m~=p && m~=1 && p~=1) || (n~=q && n~=1 && q~=1)
    error("Arrays have incompatible sizes for this operation.")
end
mC = p + (m-p)*(p==1);  nC = q + (n-q)*(q==1);
if mC==0 || nC==0
    C = quadPoly(zeros(mC,nC));
    return
end
A = repmat(A, 1+(mC-1)*(m==1), 1+(nC-1)*(n==1));
B = repmat(B, 1+(mC-1)*(p==1), 1+(nC-1)*(q==1));

% Deal with numeric .* quadPoly or quadPoly .* numeric case:
% scale each block of the coefficient matrix by the matching element
if isnumeric(A) || isnumeric(B)
    if isnumeric(A)
        [A,B] = deal(B,A);
    end
    ds = prod(cellfun(@numel,A.Zs));    dt = prod(cellfun(@numel,A.Zt));
    Cmat = kron(sparse(double(B)),ones(ds,dt)).*A.C;
    C = quadPoly(Cmat, A.Zs, A.Zt, [mC,nC], A.ns, A.nt);
    return
end

% Deal with quadPoly .* quadPoly case: compute C(:) = diag(A(:))*B(:),
% where vec(diag(A(:))) = S*A(:) for a selection matrix S
N = mC*nC;
S = sparse((0:N-1)*N+(1:N), 1:N, 1, N^2, N);
Adiag = reshape(S*reshape(A,N,1), N, N);
C = reshape(Adiag*reshape(B,N,1), mC, nC);

end
