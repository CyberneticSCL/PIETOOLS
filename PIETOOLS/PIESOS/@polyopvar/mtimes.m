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

% Make sure the inner dimensions match
if size(A,2)~=size(B,1) && ~all(size(A)==1) && ~all(size(B)==1)
    error("Incorrect dimensions for matrix multiplication.")
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
    C.C = A*B.C;
    C.dims.out = C.C.dims.out;
    C.vars.out = C.C.vars.out;
    C.dom.out = C.C.dom.out;
    C.depmat.out = C.C.depmat.out;
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