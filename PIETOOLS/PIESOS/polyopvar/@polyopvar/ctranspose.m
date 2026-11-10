function B = ctranspose(A)
% B = CTRANSPOSE(A) returns the 'polyopvar' object B = A' representing the
% L2-adjoint of the distributed polynomial A, so that the product of A'
% with a distributed polynomial V is the inner product of A and V,
%   A'*V = innerprod(A,V),
% and A'*P*V = innerprod(A,P,V) for a weight P, see 'mtimes'.
%
% INPUTS
% - A:      m x 1 'polyopvar' object, or an adjoint (A=U') with the
%           identity as weight;
%
% OUTPUTS
% - B:      'polyopvar' object storing the same polynomial as A, but with
%           B.is_adjoint = ~A.is_adjoint, and B.weight = 1;
%
% NOTES
% - An adjoint stores the polynomial U that it is the adjoint of, and not
%   the polynomial U'. Most operations are therefore not supported for an
%   adjoint, only the products U'*V and U'*P, sums of adjoints with the
%   same weight, the negative, and the horizontal concatenation of
%   adjoints, see 'mtimes', 'plus' and 'horzcat'.
% - The transpose of an object with a weight W other than the identity,
%   e.g. the transpose (U'*P)' of an adjoint, is not supported: an error
%   is thrown, whether or not A is an adjoint.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - ctranspose
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
% DJ, 11/07/2026: Initial coding;

if ~A.chkval
    A = initialize(A);
end

% Transposition of polynomial with a weight is not supported
if ~is_eye(A.weight,size(A.C,1))
    error("The transpose of an object with a weight other than the identity, e.g. (u'*P)', is not supported.")
end

B = A;
B.is_adjoint = ~A.is_adjoint;
% The weight of B is 1, also if that of A is the identity eye(n)
B.weight = 1;

end
