function C = mtimes(A,B)
% C = MTIMES(A,B) returns the 'polyopvar' object C representing the product
% A*B of two distributed polynomials, or of a constant with a distributed
% polynomial.
%
% INPUTS
% - A:     'polyopvar' object, or a constant (e.g. 'double' array) such
%          that A*B.C is defined, see @tensopmat/mtimes.m;
% - B:     'polyopvar' object, or a numeric scalar;
%
% OUTPUTS
% - C:     'polyopvar' object representing A*B;
%
% NOTES
% - The product of two 'polyopvar' objects is only supported if one of
%   them is scalar, and is computed using times.
% - Right-multiplication of a 'polyopvar' object A by a constant is only
%   supported for a numeric scalar B.
% - Left-multiplication of B by a non-scalar 'double' matrix A returns a
%   polynomial with a single row-block, each row of which is a linear
%   combination of the rows of B. If the rows of B map to different output
%   spaces (e.g. R and L2[s], or L2[s1] and L2[s2]), they are first lifted
%   to a common output space, e.g. elements of R become elements of L2[s1]
%   that are constant in s1.
% - Multiplication of a non-scalar matrix with a monomial on which an
%   operator with several factors acts (e.g. x0*u) is not supported, see
%   @tensopvar/mtimes.m.
% - If A is an adjoint (A=U'*W, see 'ctranspose'), A*B is
%   - for a 'polyopvar' B, the inner product innerprod(U,W,B), a scalar
%     'polyopvar', so that U'*B = innerprod(U,B);
%   - for a 'double', 'polynomial' or 'dpvar' B with as many rows as A has
%     columns, or a scalar, the adjoint U'*W*B, with weight W*B (B for the
%     unit weight W=1), so that U'*P*V = innerprod(U,P,V). A scalar weight
%     P represents P*eye(n).
%   A product with an adjoint as the right factor, e.g. U*V', is not
%   supported.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mtimes
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
% DJ, 10/15/2026: Initial coding;
% DJ, 11/04/2026: Left-multiplication by a matrix lifts rows with different
%                   output spaces to a common space;
% DJ, 11/06/2026: Output space of C follows from C.C (see set.C);
% DJ, 11/07/2026: Support products with an adjoint, see 'ctranspose';

% Make sure the inner dimensions match
if size(A,2)~=size(B,1) && ~all(size(A)==1) && ~all(size(B)==1)
    error("Incorrect dimensions for matrix multiplication.")
end

% Handle an adjoint separately: it is not a polynomial in the states
if isa(B,'polyopvar') && B.is_adjoint
    error("Multiplication with an adjoint as the right factor, e.g. u*v', is not supported.")
elseif isa(A,'polyopvar') && A.is_adjoint
    % The polynomial U of A=U'*W, without the weight W
    Au = A;     Au.is_adjoint = false;      Au.weight = 1;
    if isa(B,'polyopvar')
        % <U,W*B>: innerprod checks the dimensions of W, as the check above
        % skips a scalar A
        C = innerprod(Au,A.weight,B);
    elseif (isnumeric(B) || isa(B,'polynomial') || isa(B,'dpvar')) && (all(size(B)==1) || size(B,1)==size(A,2))
        % U'*W*B stores the weight W*B, for the inner product with a later
        % factor. Check the dimensions of B here, as the check above skips
        % a scalar A
        C = A;
        C.weight = A.weight*B;
    else
        error("An adjoint can only be multiplied with a 'polyopvar' object, or with a scalar or matrix of type 'double', 'polynomial' or 'dpvar' with as many rows as the adjoint has columns.")
    end
    return
end

% Distinguish three cases
if isa(A,'polyopvar') && isa(B,'polyopvar')
    % Support multiplication of polynomials only if one is scalar for now
    if all(size(A)==1) || all(size(B)==1)
        C = A.*B;
        return
    else
        error("Multiplication of matrix-valued distributed polynomials is currently not supported.")
    end
elseif ~isa(A,'polyopvar')
    % Multiplication of constant with polynomial
    % A row of A*B combines the rows of B, so these must map to the same
    % output space: lift them to a common one if they differ, as in plus.
    % An operator A acts on the output space of each row as it is: no lift.
    if isa(A,'double') && ~all(size(A)==1)
        if ~B.chkval
            B = initialize(B);
        end
        if any(any(diff(B.depmat.out,1,1)))                                 % DJ, 11/04/2026
            B = lift_output(B,B.vars,B.dom,[],max(B.depmat.out,[],1));
        end
    end
    C = B;
    % The output space of C is set by the assignment of C.C, see set.C
    C.C = A*B.C;
    return
elseif ~isa(B,'polyopvar')
    % Multiplication of polynomial with constant
    if isnumeric(B) && all(size(B)==1)
        C = A;
        C.C = A.C*B;
        return
    else
        error("Right-multiplication of polynomials is currently not supported.")
    end
end

end