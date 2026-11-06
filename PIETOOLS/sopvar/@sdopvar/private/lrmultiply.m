function Cout = lrmultiply(L, C, R)
% COUT = LRMULTIPLY(L,C,R) left- and right-multiplies a vectorized
% decision-variable-dependent coefficient by constant matrices, computing
%
%   Cout = vec(L*C*R),    C = A + B'*d,
%
% without forming kron(R',L) explicitly for the constant part.
%
% INPUTS
% - L:  sparse matrix, left factor;
% - C:  struct with fields 'A' (vec of the constant part) and 'B' (rows
%       indexed by decision variables, columns by vec of the matrix);
% - R:  sparse matrix, right factor;
%
% OUTPUTS
% - Cout:   struct with fields 'A' and 'B', representing vec(L*C*R);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lrmultiply
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
% Initial coding SS; vectorization AT, 09/08/2026
% MMP, 11/06/2026: Comment pass (maintainer, 11/06/2026): commented-out
%                  earlier implementations removed, step comments; the
%                  preallocation of Bout before the product is dead and is
%                  left in place; code unchanged.


A = C.A;
B = C.B;

% Constant part: vec(L*X*R) = (R' kron L)*vec(X) with X = reshape(A),
% without forming the Kronecker product
X = reshape(A,size(L,2),size(R,1));
Y = L*X*R;  % (R^T\otimes L)*A = vec(L*X*R), X = reshape(A)
Aout = Y(:);

Bout = sparse(size(B,1),size(L,1)*size(R,2));

% Decision part, one row vec(X_i)' of B per variable: the rows map by the
% transpose, B*(R kron L')
LRkron = kron(R, L');
Bout = B*LRkron;

Cout = struct('A',Aout,'B', Bout);
end