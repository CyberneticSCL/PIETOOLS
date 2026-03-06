function C = mtimes(A,B)
% C = mtimes(A,B) returns the 'polyopvar' object C representing the product
% of two polyopvar objects.
%
% INPUTS
% - A:     polyopvar object; 
% - B:     polyopvar object;
%
% OUTPUTS
% - C:     polyopvar object representing A*B.
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
% CR, 01/28/2026: Initial coding
% DJ, 03/03/2026: Update using 'otimes' function, get rid of duplicate
%                   monomials
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make sure the inner dimensions match
if size(A,2)~=size(B,2) && ~all(size(A)==1) && ~all(size(B)==1)
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
    return
elseif ~isa(B,'polyopvar')
    % Multiplication of polynomial with constant
    C = A;
    C.C = A.C*B;
    return
end

end