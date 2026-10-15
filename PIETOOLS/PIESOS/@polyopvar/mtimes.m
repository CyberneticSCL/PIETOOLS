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